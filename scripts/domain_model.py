"""
DomainModel class for bone healing simulation using element-based geometry.
"""

import numpy as np
from pathlib import Path
from mesa import Model
from mesa.time import RandomActivation

from endothelial_cell_agent import EndothelialCellAgent

from element_agent_optimized import ElementAgent
from mesh import MeshGeometry
from reaction import react_element, decay_unoccupied, production_rates
from bone_healing_model_optimized import hypoxia_regulation
from parameters import DT_HOURS


class DomainModel(Model):
    def __init__(self, nodes, elements, params, seed=None):
        """
        Initialize the DomainModel using an element-based geometry.
        :param nodes: Dictionary of nodes with their coordinates.
        :param elements: Dictionary of elements with their node IDs.
        :param params: Model parameters for the bone healing equations.
        :param seed: Seed for the model's random number generator (self.random).
                     All stochastic choices draw from self.random, so a given seed
                     reproduces a run exactly. Mesa reads it in Model.__new__.
        """
        super().__init__()
        # Parcels are activated in a fresh random order every step (seeded through self.random)
        self.schedule = RandomActivation(self)
        self.nodes = nodes
        self.elements = elements
        self.params = params
        self.time = 0

        self.vessel_segments = []
        self.element_agents = {element_id: [] for element_id in self.elements.keys()}
        self.oxygen_field = {element_id: 0.0 for element_id in self.elements}
        self.vessel_loops = {}
        self.vessel_element_ids = set()

        self.enable_EC = False
        self.debug = False

        # Compute centroids for elements
        self.element_centroids = {
            element_id: tuple(
                np.mean([nodes[node_id] for node_id in node_ids], axis=0)
            )
            for element_id, node_ids in elements.items()
        }
        self.coord_to_element = {(int(x), int(y)): eid for eid, (x, y) in self.element_centroids.items()}

        # Mesh geometry: element areas, edge adjacency, finite-volume transport operator
        self.geometry = MeshGeometry(nodes, elements)
        self.area_ratio = dict(zip(self.geometry.ids, self.geometry.area_ratio))
        self.max_parcels = int(self.params.get("max_parcels_per_element", 3))
        self.n_clipped = 0

        # Initialize debris field
        self.debris_field = {}

        # NEW: coefficient from params (default = 1.0)
        debris_coeff = float(self.params.get("debris_coeff", 1.0))

        for element_id, centroid in self.element_centroids.items():
            centroid_x, centroid_y = centroid

            if centroid_x < 2.0:
                self.debris_field[element_id] = 200 * debris_coeff
            elif 2.0 <= centroid_x <= 7.5:
                decay_constant = 1.0
                self.debris_field[element_id] = 200 * debris_coeff * np.exp(-decay_constant * (centroid_x - 2.0))
            if 2.5 <= centroid_y <= 10 or -2.5 >= centroid_y >= -10:
                max_y = 10
                min_y = 2.5
                self.debris_field[element_id] = self.debris_field[element_id] * ((max_y - abs(centroid_y)) / (max_y - min_y))
                self.debris_field[element_id] = max(0, self.debris_field[element_id])
            else:
                self.debris_field[element_id] = self.debris_field[element_id]

        # Initialize cytokine fields
        self.cytokine_fields = {
            "c1": {},
            "c2": {},
            "c3": {},
            "c4": {},
        }

        radii = {
            "c1": 3.5,
            "c2": 2.0,
            "c3": 1.5,
            "c4": 4.0,
        }

        for element_id, node_ids in elements.items():
            centroid_x = np.mean([nodes[node_id][0] for node_id in node_ids])
            centroid_y = np.mean([nodes[node_id][1] for node_id in node_ids])
            distance_from_center = np.sqrt(centroid_x**2 + centroid_y**2)

            for cytokine, radius in radii.items():
                if distance_from_center <= radius:
                    self.cytokine_fields[cytokine][element_id] = params.get(f"init_{cytokine}", {
                        "c1": 0.050157, "c2": 0.0, "c3": 0.0, "c4": 0.02
                    }[cytokine])
                else:
                    self.cytokine_fields[cytokine][element_id] = 0.0

        self._initialize_agents()

        # Neighbours are elements sharing an edge (not merely a corner node)
        self.neighbor_cache = self.geometry.neighbors

    def snapshot_outputs(self):
        """Capture current state for validation or analysis."""
        counts = {
            "PMN": 0.0,
            "M0": 0.0,
            "M1": 0.0,
            "M2": 0.0,
            "MSC": 0.0,
            "EC": 0
        }

        for a in self.schedule.agents:
            if isinstance(a, ElementAgent):
                counts["PMN"] += a.state[0]
                counts["M0"] += a.state[1]
                counts["M1"] += a.state[2]
                counts["M2"] += a.state[3]
                counts["MSC"] += a.state[7]
            elif isinstance(a, EndothelialCellAgent):
                counts["EC"] += 1

        totals = self.cytokine_amounts()
        return {"counts": counts, "cytokines": totals}

    def is_closed_loop(self, path, tol=0.1):
        """Check if a vessel path forms a closed loop."""
        if len(path) < 3:
            return False
        x0, y0 = path[0]
        x1, y1 = path[-1]
        return np.linalg.norm(np.array([x0, y0]) - np.array([x1, y1])) < tol

    def remove_parcel(self, agent):
        """Retire a parcel: remove it from the scheduler and from its element's occupancy list."""
        if agent in self.element_agents[agent.element_id]:
            self.element_agents[agent.element_id].remove(agent)
        self.schedule.remove(agent)

    def move_agent_to(self, agent, target_element):
        """Move an agent to a new element and update references."""
        if agent in self.element_agents[agent.element_id]:
            self.element_agents[agent.element_id].remove(agent)

        agent.element_id = target_element
        agent.centroid = self.element_centroids[target_element]
        self.element_agents[target_element].append(agent)

    def _initialize_agents(self):
        """Initialize agents based on element centroids."""

        # NEW: allow JSON-controlled initial agent counts with defaults
        num_pmn_agents = int(self.params.get("init_PMN_agents", 100))
        num_m0_agents  = int(self.params.get("init_M0_agents", 1))
        num_cm_agents  = int(self.params.get("init_MSC_agents", 2))

        # Add PMN agents
        for i in range(num_pmn_agents):
            eligible_elements = [
                (element_id, centroid) for element_id, centroid in self.element_centroids.items()
                if -2.5 <= centroid[0] <= 6.0
            ]

            weights = []
            for element_id, centroid in eligible_elements:
                if -2.5 <= centroid[1] <= 2.5:
                    weights.append(3)
                else:
                    weights.append(1)

            chosen_element = self.random.choices(eligible_elements, weights=weights, k=1)[0]
            element_id, centroid = chosen_element

            initial_conditions = [
                10, 0, 0, 0,
                self.cytokine_fields["c1"][element_id],
                self.cytokine_fields["c2"][element_id],
                self.cytokine_fields["c3"][element_id], 0,
                self.cytokine_fields["c4"][element_id]
            ]
            agent = ElementAgent(f"PMN-{i}", self, initial_conditions, self.params, element_id, centroid)
            self.schedule.add(agent)
            self.element_agents[element_id].append(agent)

        # Add M0 agents
        for j in range(num_m0_agents):
            eligible_elements = [
                element_id for element_id, centroid in self.element_centroids.items()
                if -2.0 <= centroid[0] <= 4
            ]
            element_id = self.random.choice(eligible_elements)
            centroid = self.element_centroids[element_id]
            initial_conditions = [
                0,
                self.params.get("init_M0", 1),
                self.params.get("init_M1", 0),
                self.params.get("init_M2", 0),
                self.cytokine_fields["c1"][element_id],
                self.cytokine_fields["c2"][element_id],
                self.cytokine_fields["c3"][element_id],
                0,
                self.cytokine_fields["c4"][element_id],
            ]
            agent = ElementAgent(f"M0-{j}", self, initial_conditions, self.params, element_id, centroid)
            self.schedule.add(agent)
            self.element_agents[element_id].append(agent)

        # Add MSC agents
        msc_counter = 0
        # Layout 0 (v1.0.1): two parcels of 5 cell-equivalents at the outer margin of the gap.
        # Layout 1: the same total progenitor amount distributed over parcels in the periosteal
        # cambium layer along the trans-cortex surface (first 0.3 mm of callus, |y| > 1 mm).
        # Layout 2: half of the total in the periosteal cambium layer, half in the two gap parcels
        # (periosteal and marrow-side progenitor sources); total amount unchanged.
        layout = int(self.params.get("init_MSC_layout", 0))
        gap_amount = 5.0 if layout != 2 else 2.5
        if layout in (1, 2):
            total_msc = 5.0 * num_cm_agents * (1.0 if layout == 1 else 0.5)
            n_parcels = int(self.params.get("init_MSC_periosteal_parcels", 10))
            cambium = [eid for eid, (x, y) in self.element_centroids.items()
                       if 1.88 <= x <= 2.18 and 1.0 < abs(y) <= 10.0]
            for element_id in self.random.sample(cambium, n_parcels):
                c = self.cytokine_fields
                initial_conditions = [0, 0, 0, 0, c["c1"][element_id], c["c2"][element_id],
                                      c["c3"][element_id], total_msc / n_parcels, c["c4"][element_id]]
                agent = ElementAgent(f"MSC-{msc_counter}", self, initial_conditions, self.params,
                                     element_id, self.element_centroids[element_id])
                self.schedule.add(agent)
                self.element_agents[element_id].append(agent)
                msc_counter += 1
            if layout == 1:
                num_cm_agents = 0
        half_num_cm_agents = num_cm_agents // 2

        # Group 1: MSCs in range (-3, 0)
        for k in range(half_num_cm_agents):
            eligible_elements = [
                element_id for element_id, centroid in self.element_centroids.items()
                if 1.2 <= centroid[0] <= 1.6 and -3 <= centroid[1] < 0
            ]

            element_id = self.random.choice(eligible_elements)
            centroid = self.element_centroids[element_id]
            initial_conditions = [
                0, 0, 0, 0,
                self.cytokine_fields["c1"][element_id],
                self.cytokine_fields["c2"][element_id],
                self.cytokine_fields["c3"][element_id],
                gap_amount,
                self.cytokine_fields["c4"][element_id],
            ]
            agent = ElementAgent(f"MSC-{msc_counter}", self, initial_conditions, self.params, element_id, centroid)
            self.schedule.add(agent)
            self.element_agents[element_id].append(agent)
            msc_counter += 1

        # Add initial EC agents if enabled
        if self.enable_EC:
            ec_counter = 0
            eligible_elements = []

            for element_id, centroid in self.element_centroids.items():
                x, y = centroid
                if (-2 < x < 1.4 and (-1 <= y <= -0.5 or 0.5 <= y <= 1)) or (1.5 < x < 1.7 and -2 < y < 2):
                    eligible_elements.append((element_id, centroid))

            num_to_select = int(0.7 * len(eligible_elements))
            selected_elements = self.random.sample(eligible_elements, num_to_select)

            for element_id, centroid in selected_elements:
                ec_agent = EndothelialCellAgent(f"EC-{ec_counter}", self, element_id, centroid)
                self.schedule.add(ec_agent)
                self.element_agents[element_id].append(ec_agent)
                ec_counter += 1

        # Group 2: MSCs in range (0, 3)
        for k in range(num_cm_agents - half_num_cm_agents):
            eligible_elements = [
                element_id for element_id, centroid in self.element_centroids.items()
                if 1.2 <= centroid[0] <= 1.6 and 0 <= centroid[1] <= 3
            ]

            element_id = self.random.choice(eligible_elements)
            centroid = self.element_centroids[element_id]
            initial_conditions = [
                0, 0, 0, 0,
                self.cytokine_fields["c1"][element_id],
                self.cytokine_fields["c2"][element_id],
                self.cytokine_fields["c3"][element_id],
                gap_amount,
                self.cytokine_fields["c4"][element_id],
            ]
            agent = ElementAgent(f"MSC-{msc_counter}", self, initial_conditions, self.params, element_id, centroid)
            self.schedule.add(agent)
            self.element_agents[element_id].append(agent)
            msc_counter += 1

    def find_branch_location(self, element_id):
        """Find a neighboring element with no ECs for branching."""
        neighbors = self.get_neighbors(element_id)
        viable = [
            eid for eid in neighbors
            if all(not isinstance(a, EndothelialCellAgent) for a in self.element_agents[eid])
        ]
        return self.random.choice(viable) if viable else None

    def link_sprouts(self, ec1, ec2):
        """Store loop information between two ECs."""
        key = (ec1.unique_id, ec2.unique_id)
        self.vessel_loops[key] = ec1.sprout_path + ec2.sprout_path
        combined_path = ec1.sprout_path + ec2.sprout_path
        new_segment = {
            "path": combined_path,
            "maturity": "mature",
            "age": 0,
            "last_high_VEGF_time": self.time
        }
        self.vessel_segments.append(new_segment)

    def step_agents(self):
        """Step all agents."""
        def agent_step(agent):
            agent.step()

    def find_available_element(self, element_id, visited=None, depth_limit=50):
        """Find an available element for agent placement within a depth limit."""
        if visited is None:
            visited = set()

        visited.add(element_id)

        if len(self.element_agents[element_id]) < self.max_parcels:
            return element_id

        if depth_limit <= 0:
            return None

        neighbors = self.get_neighbors(element_id)
        for neighbor in neighbors:
            if neighbor not in visited:
                available = self.find_available_element(neighbor, visited, depth_limit - 1)
                if available is not None:
                    return available

        return None

    def spawn_agent(self, element_id, initial_conditions):
        """Spawn a new agent in a neighboring element."""
        target_element_id = self.find_available_element(element_id)

        if target_element_id is None:
            return

        if initial_conditions[1] > 0:
            agent_type = "M0"
        elif initial_conditions[2] > 0:
            agent_type = "M1"
        elif initial_conditions[3] > 0:
            agent_type = "M2"
        elif initial_conditions[7] > 0:
            agent_type = "MSC"
        else:
            agent_type = "Unknown"

        agent_count = sum(1 for agent in self.schedule.agents if agent_type in agent.unique_id)
        new_id = f"{agent_type}-{agent_count + 1}"

        centroid = self.element_centroids[target_element_id]
        new_agent = ElementAgent(new_id, self, initial_conditions, self.params, target_element_id, centroid)
        self.schedule.add(new_agent)
        self.element_agents[target_element_id].append(new_agent)

    def update_oxygen_field(self):
        """Update oxygen levels based on EC delivery, diffusion, and consumption."""
        new_oxygen = self.oxygen_field.copy()

        for segment in self.vessel_segments:
            if segment["maturity"] != "mature":
                continue
            for x, y in segment["path"]:
                coord = (int(x), int(y))
                element_id = self.coord_to_element.get(coord)
                if element_id is not None:
                    new_oxygen[element_id] += 0.2

        diffused_oxygen = new_oxygen.copy()
        for element_id, value in new_oxygen.items():
            neighbors = self.get_neighbors(element_id)
            centroid = self.element_centroids[element_id]
            for neighbor in neighbors:
                neighbor_value = new_oxygen[neighbor]
                neighbor_centroid = self.element_centroids[neighbor]
                distance = np.linalg.norm(np.array(centroid) - np.array(neighbor_centroid))
                decay = np.exp(-distance)
                flux = self.params.get("diff_O2", 1e-2) * decay * (value - neighbor_value)
                if flux > 0:
                    diffused_oxygen[element_id] -= flux
                    diffused_oxygen[neighbor] += flux

        for agent in self.schedule.agents:
            if isinstance(agent, ElementAgent):
                eid = agent.element_id
                consumption = 0.05 * (agent.state[0] + agent.state[1] + agent.state[2] + agent.state[3] + agent.state[7])
                diffused_oxygen[eid] = max(0, diffused_oxygen[eid] - consumption)

        self.oxygen_field = diffused_oxygen

    # ------------------------------------------------------------------
    # Continuous updates (Box 1, steps 1-3)
    # ------------------------------------------------------------------
    CYTOKINES = ("c1", "c2", "c3", "c4")
    CELL_INDEX = [0, 1, 2, 3, 7]        # PMN, M0, M1, M2, MSC in ElementAgent.state
    FIELD_INDEX = [4, 5, 6, 8]          # c1, c2, c3, c4 mirrored in ElementAgent.state

    def cytokine_amounts(self):
        """Domain amount of each cytokine, sum_e a_e * c_e (in mean-element units); matches bulk dPCR."""
        return {c: float(sum(self.area_ratio[e] * v for e, v in self.cytokine_fields[c].items()))
                for c in self.CYTOKINES}

    def cytokine_production(self):
        """Domain secretion rate of each cytokine (amount per hour), summed over occupied elements."""
        by_element = {}
        for a in self.schedule.agents:
            if isinstance(a, ElementAgent):
                by_element.setdefault(a.element_id, []).append(a)
        total = np.zeros(4)
        for eid, parcels in by_element.items():
            fields = np.array([self.debris_field[eid]] + [self.cytokine_fields[k][eid] for k in self.CYTOKINES])
            cells = np.array([a.state[self.CELL_INDEX] for a in parcels], dtype=float)
            total += production_rates(fields, cells, self.params, self.area_ratio[eid],
                                      hypoxia_regulation(self.oxygen_field[eid]))
        return dict(zip(self.CYTOKINES, total.tolist()))

    def react(self, dt=DT_HOURS):
        """Element-local reaction step: parcels + cytokines + debris of each element integrated together."""
        by_element = {}
        for a in self.schedule.agents:
            if isinstance(a, ElementAgent):
                by_element.setdefault(a.element_id, []).append(a)

        for eid in self.elements:
            c = np.array([self.cytokine_fields[k][eid] for k in self.CYTOKINES])
            parcels = by_element.get(eid)
            if not parcels:
                new_c = decay_unoccupied(c, self.params, dt)
                for k, v in zip(self.CYTOKINES, new_c):
                    self.cytokine_fields[k][eid] = float(v)
                continue
            fields = np.concatenate([[self.debris_field[eid]], c])
            cells = np.array([a.state[self.CELL_INDEX] for a in parcels], dtype=float)
            new_fields, new_cells, clipped = react_element(
                fields, cells, self.params, self.area_ratio[eid],
                EC=self.element_EC_counts.get(eid, 0), PO2=self.oxygen_field[eid], dt=dt)
            self.n_clipped += clipped
            self.debris_field[eid] = float(new_fields[0])
            for k, v in zip(self.CYTOKINES, new_fields[1:]):
                self.cytokine_fields[k][eid] = float(v)
            for a, row in zip(parcels, new_cells):
                a.state = a.state.astype(float)
                a.state[self.CELL_INDEX] = row

    def transport(self, dt=DT_HOURS):
        """Finite-volume diffusion of each cytokine over the edge graph (mass-conserving, zero-flux)."""
        ids = self.geometry.ids
        for k in self.CYTOKINES:
            arr = np.array([self.cytokine_fields[k][e] for e in ids])
            arr = self.geometry.diffuse(arr, float(self.params[f"diff_{k}"]), dt)
            self.cytokine_fields[k] = dict(zip(ids, arr.tolist()))

    def mirror_fields_to_parcels(self):
        """Keep the cytokine slots of each parcel's state equal to its element's field (read-only copy)."""
        for a in self.schedule.agents:
            if isinstance(a, ElementAgent):
                for idx, k in zip(self.FIELD_INDEX, self.CYTOKINES):
                    a.state[idx] = self.cytokine_fields[k][a.element_id]

    # ------------------------------------------------------------------
    # Chemotaxis (Box 1, step 4; Methods 5.5)
    # ------------------------------------------------------------------
    def migrate_parcels(self):
        """
        Composition-weighted chemotaxis. Each parcel weights the debris, TNF-a, IL-10, TGF-b and VEGF
        fields by its own composition (PMN/M0 -> debris, M1 -> TNF-a, M2 -> IL-10, MSC -> TGF-b and VEGF),
        scores the current element and every edge neighbour with free capacity by the weighted change in
        max-normalised fields, and moves to a candidate drawn from a softmax with temperature tau.
        """
        beta_N = float(self.params.get("chemotaxis_beta_N", 0.6))
        beta_0 = float(self.params.get("chemotaxis_beta_0", 1.0))
        tau = float(self.params.get("chemotaxis_tau", 0.05))

        def normalised(field):
            m = max(field.values())
            return {e: (v / m if m > 0 else 0.0) for e, v in field.items()}

        phi = {"D": normalised(self.debris_field)}
        for k in self.CYTOKINES:
            phi[k] = normalised(self.cytokine_fields[k])

        parcels = [a for a in self.schedule.agents if isinstance(a, ElementAgent)]
        self.random.shuffle(parcels)
        for a in parcels:
            N, M0, M1, M2, C = (a.state[i] for i in self.CELL_INDEX)
            Y = N + M0 + M1 + M2 + C
            if Y <= 0:
                continue
            w = {"D": (beta_N * N + beta_0 * M0) / Y, "c1": M1 / Y, "c2": M2 / Y, "c3": C / Y, "c4": C / Y}
            here = a.element_id
            candidates = [here] + [j for j in self.get_neighbors(here)
                                   if len(self.element_agents[j]) < self.max_parcels]
            scores = np.array([sum(w[s] * (phi[s][j] - phi[s][here]) for s in w) for j in candidates])
            weights = np.exp((scores - scores.max()) / tau)
            target = self.random.choices(candidates, weights=weights.tolist(), k=1)[0]
            if target != here:
                self.move_agent_to(a, target)

    def get_neighbors(self, element_id):
        """
        Get neighboring elements based on shared nodes (using cached lookup).

        This is O(1) lookup instead of O(N) iteration through all elements.

        :param element_id: The ID of the current element
        :return: List of neighboring element IDs
        """
        return self.neighbor_cache.get(element_id, [])

    def get_element_id_from_centroid(self, target_centroid):
        """Find the closest element whose centroid matches the given (x, y)."""
        for eid, centroid in self.element_centroids.items():
            if np.allclose(centroid, target_centroid, atol=1e-5):
                return eid
        raise ValueError(f"No matching element found for centroid {target_centroid}")

    def step(self):
        """One global step of DT_HOURS (Methods Box 1)."""
        if self.enable_EC:
            self.element_EC_counts = {
                eid: sum(1 for a in agents if isinstance(a, EndothelialCellAgent))
                for eid, agents in self.element_agents.items()
            }
        else:
            self.element_EC_counts = {eid: 0 for eid in self.elements.keys()}
        self.update_oxygen_field()

        # 2-3. reaction and transport (Lie splitting). n_split > 1 alternates them n times per hour
        # with dt = 1/n h; used only for the splitting-error test (Methods 5.6.3). Default 1.
        n_split = int(self.params.get("n_split", 1))
        for _ in range(n_split):
            self.react(DT_HOURS / n_split)      # element-local reaction (cells, cytokines, debris)
            self.transport(DT_HOURS / n_split)  # cytokine transport
        self.mirror_fields_to_parcels()
        self.population_totals = np.sum(
            [a.state for a in self.schedule.agents if isinstance(a, ElementAgent)], axis=0)
        self.schedule.step()              # 4a. parcel events in random order: budding (ElementAgent.step)
        self.migrate_parcels()            # 4b. chemotaxis in random order

        VEGF_THRESHOLD = 0.01
        PRUNING_DELAY = 6

        for segment in self.vessel_segments:
            segment["age"] += 1
            try:
                avg_vegf = np.mean([self.cytokine_fields["c4"][int(y)][int(x)] for x, y in segment["path"]])
            except:
                continue

            if avg_vegf > VEGF_THRESHOLD:
                segment["last_high_VEGF_time"] = self.time
            else:
                time_since_high = self.time - segment["last_high_VEGF_time"]
                if time_since_high > PRUNING_DELAY and segment["maturity"] == "immature":
                    segment["maturity"] = "pruned"
        self.time += 1