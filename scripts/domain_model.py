"""
DomainModel class for bone healing simulation using element-based geometry.
"""

import numpy as np
import random
from pathlib import Path
from mesa import Model
from mesa.time import SimultaneousActivation

#from element_agent import ElementAgent
from endothelial_cell_agent import EndothelialCellAgent
#from bone_healing_model import full_bone_healing_model

from element_agent_optimized import ElementAgent
from bone_healing_model_optimized import full_bone_healing_model


class DomainModel(Model):
    def __init__(self, nodes, elements, params):
        """
        Initialize the DomainModel using an element-based geometry.
        :param nodes: Dictionary of nodes with their coordinates.
        :param elements: Dictionary of elements with their node IDs.
        :param params: Model parameters for the bone healing equations.
        """
        super().__init__()
        self.schedule = SimultaneousActivation(self)
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

        # Pre-calculate neighbor relationships for O(1) lookup instead of O(N) iteration
        self.neighbor_cache = self._build_neighbor_cache()

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

        totals = {"c1": 0.0, "c2": 0.0, "c3": 0.0, "c4": 0.0}

        for a in self.schedule.agents:
            if isinstance(a, ElementAgent):
                counts["PMN"] += a.state[0]
                counts["M0"] += a.state[1]
                counts["M1"] += a.state[2]
                counts["M2"] += a.state[3]
                counts["MSC"] += a.state[7]

                totals["c1"] += a.state[4]
                totals["c2"] += a.state[5]
                totals["c3"] += a.state[6]
                totals["c4"] += a.state[8]

            elif isinstance(a, EndothelialCellAgent):
                counts["EC"] += 1

        return {"counts": counts, "cytokines": totals}

    def _build_neighbor_cache(self):
        """
        Pre-calculate all neighbor relationships once during initialization.
        This is much faster than computing neighbors on-the-fly every time.

        :return: Dictionary mapping element_id -> list of neighbor element_ids
        """
        neighbor_cache = {}

        for element_id, node_ids in self.elements.items():
            neighbors = set()
            current_nodes = set(node_ids)

            for other_element, other_nodes in self.elements.items():
                if element_id != other_element and current_nodes.intersection(other_nodes):
                    neighbors.add(other_element)

            neighbor_cache[element_id] = list(neighbors)

        return neighbor_cache

    def is_closed_loop(self, path, tol=0.1):
        """Check if a vessel path forms a closed loop."""
        if len(path) < 3:
            return False
        x0, y0 = path[0]
        x1, y1 = path[-1]
        return np.linalg.norm(np.array([x0, y0]) - np.array([x1, y1])) < tol

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
        half_num_cm_agents = num_cm_agents // 2
        msc_counter = 0

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
                5,
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
            selected_elements = random.sample(eligible_elements, num_to_select)

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
                5,
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

        if len(self.element_agents[element_id]) < 3:
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

    def update_cytokine_fields(self):
        """Update the cytokine fields based on the agents' states."""
        for cytokine in self.cytokine_fields:
            for element_id in self.cytokine_fields[cytokine]:
                self.cytokine_fields[cytokine][element_id] = 0.0

        for agent in self.schedule.agents:
            if not isinstance(agent, ElementAgent):
                continue
            element_id = agent.element_id
            self.cytokine_fields["c1"][element_id] += agent.state[4]
            self.cytokine_fields["c2"][element_id] += agent.state[5]
            self.cytokine_fields["c3"][element_id] += agent.state[6]
            self.cytokine_fields["c4"][element_id] += agent.state[8]

    def update_debris_field(self):
        """Update the debris field based on agent interactions."""
        updated_debris_field = self.debris_field.copy()

        for agent in self.schedule.agents:
            if not isinstance(agent, ElementAgent):
                continue
            element_id = agent.element_id
            if updated_debris_field[element_id] > 0:
                D = updated_debris_field[element_id]
                M0, M1, M2, PMN = agent.state[1], agent.state[2], agent.state[3], agent.state[0]
                R_D = D / (self.params["a_ed"] + D)
                debris_consumed = R_D * (5 * self.params["k_e0"] * M0 +
                                         3 * self.params["k_e1"] * M1 +
                                         2 * self.params["k_e2"] * M2 +
                                         10 * self.params["k_e_pmn"] * PMN)
                updated_debris_field[element_id] = max(0, updated_debris_field[element_id] - debris_consumed)

        self.debris_field = updated_debris_field

    def diffuse_cytokines(self):
        """Diffuse cytokines across the element-based domain with exponential decay."""
        updated_fields = {cytokine: {element_id: 0 for element_id in self.elements.keys()}
                          for cytokine in ["c1", "c2", "c3", "c4"]}

        for agent in self.schedule.agents:
            if not isinstance(agent, ElementAgent):
                continue
            element_id = agent.element_id
            updated_fields["c1"][element_id] += agent.state[4]
            updated_fields["c2"][element_id] += agent.state[5]
            updated_fields["c3"][element_id] += agent.state[6]
            self.cytokine_fields["c4"][element_id] += agent.state[8]

        for cytokine, field in self.cytokine_fields.items():
            diff_field = field.copy()

            for element_id, value in field.items():
                neighbors = self.get_neighbors(element_id)
                centroid = self.element_centroids[element_id]

                for neighbor in neighbors:
                    neighbor_centroid = self.element_centroids[neighbor]
                    distance = np.linalg.norm(np.array(centroid) - np.array(neighbor_centroid))
                    decay_factor = np.exp(-distance)

                    diffusion = self.params[f"diff_{cytokine}"] * decay_factor * (value - field[neighbor])
                    if diffusion > 0:
                        diff_field[neighbor] += diffusion
                        diff_field[element_id] -= diffusion

            diff_field = {element_id: max(0, concentration) for element_id, concentration in diff_field.items()}
            updated_fields[cytokine] = diff_field

        self.cytokine_fields = updated_fields

        for agent in self.schedule.agents:
            element_id = agent.element_id
            agent.state[4] = self.cytokine_fields["c1"][element_id]
            agent.state[5] = self.cytokine_fields["c2"][element_id]
            agent.state[6] = self.cytokine_fields["c3"][element_id]
            agent.state[8] = self.cytokine_fields["c4"][element_id]

    def perform_migration(self):
        """Perform migration based on cytokine and debris gradients."""
        gradients = {cytokine: {} for cytokine in ["debris", "c1", "c2", "c3", "c4"]}
        for element_id, value in self.debris_field.items():
            neighbors = self.get_neighbors(element_id)
            gradients["debris"][element_id] = sum(self.debris_field[neighbor] - value for neighbor in neighbors)

        for cytokine in ["c1", "c2", "c3", "c4"]:
            for element_id, value in self.cytokine_fields[cytokine].items():
                neighbors = self.get_neighbors(element_id)
                gradients[cytokine][element_id] = sum(self.cytokine_fields[cytokine][neighbor] - value
                                                      for neighbor in neighbors)

        for agent in self.schedule.agents:
            if not isinstance(agent, ElementAgent):
                continue
            curr_element_id = agent.element_id

            weighted_gradients = {
                "debris": agent.state[0] * 0.6 + agent.state[1] * 1.0,
                "c1": agent.state[2] * 1.0,
                "c2": agent.state[3] * 1.0,
                "c3": agent.state[7] * 1.0,
                "c4": agent.state[7] * 1.0,
            }

            max_gradient_element = None
            max_gradient_value = float("-inf")

            for neighbor in self.get_neighbors(curr_element_id):
                gradient_value = (
                    weighted_gradients["debris"] * gradients["debris"].get(neighbor, 0) +
                    weighted_gradients["c1"] * gradients["c1"].get(neighbor, 0) +
                    weighted_gradients["c2"] * gradients["c2"].get(neighbor, 0) +
                    weighted_gradients["c3"] * gradients["c3"].get(neighbor, 0) +
                    weighted_gradients["c4"] * gradients["c4"].get(neighbor, 0)
                )

                if gradient_value > max_gradient_value:
                    max_gradient_value = gradient_value
                    max_gradient_element = neighbor

            if max_gradient_element is not None and max_gradient_element != curr_element_id:
                agent.migrate(max_gradient_element)

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
        """Perform one simulation step."""
        if self.time % 1 == 0:
            self.update_debris_field()
            self.update_cytokine_fields()
            self.update_oxygen_field()

                    # NEW: Pre-calculate EC counts per element (saves recalculating for each agent)
            if self.enable_EC:
                self.element_EC_counts = {
                    eid: sum(1 for a in agents if isinstance(a, EndothelialCellAgent))
                    for eid, agents in self.element_agents.items()
                }
            else:
                self.element_EC_counts = {eid: 0 for eid in self.elements.keys()}

            # DEBUG LOOP HERE (if debug and time % 8)

            self.perform_migration()
            self.schedule.step()



            if self.debug and self.time % 8 == 0:
                max_polarization_M0_M1 = 0
                max_polarization_M0_M2 = 0
                max_polarization_M1_M2 = 0
                max_dM0_dt = 0
                max_dM1_dt = 0
                max_dM2_dt = 0
                max_dD_dt = 0

                for agent in self.schedule.agents:
                    if not isinstance(agent, ElementAgent):
                        continue
                    element_id = agent.element_id

                    if self.enable_EC:
                        EC = sum(1 for a in self.element_agents[element_id] if isinstance(a, EndothelialCellAgent))
                    else:
                        EC = 0

                    D = self.debris_field[element_id]
                    PO2 = self.oxygen_field[element_id]
                    variables = agent.state

                    dPMN_dt, dM0_dt, dM1_dt, dM2_dt, _, _, _, _, _ = full_bone_healing_model(0, variables, self.params, D, EC, PO2)
                    P_M0_to_M1 = self.params["k01"] * (variables[5] / (self.params["a01"] + variables[5]))
                    P_M0_to_M2 = self.params["k02"] * (variables[6] / (self.params["a02"] + variables[6]))
                    P_M1_to_M2 = self.params["k12"] * (variables[6] / (self.params["a_M1_to_M2"] + variables[6]))

                    max_polarization_M0_M1 = max(max_polarization_M0_M1, P_M0_to_M1)
                    max_polarization_M0_M2 = max(max_polarization_M0_M2, P_M0_to_M2)
                    max_polarization_M1_M2 = max(max_polarization_M1_M2, P_M1_to_M2)
                    max_dM0_dt = max(max_dM0_dt, dM0_dt)
                    max_dM1_dt = max(max_dM1_dt, dM1_dt)
                    max_dM2_dt = max(max_dM2_dt, dM2_dt)

                print(
                    f"Iteration {self.time}: Max P_M0->M1: {max_polarization_M0_M1:.5f}, "
                    f"Max P_M0->M2: {max_polarization_M0_M2:.5f}, Max P_M1->M2: {max_polarization_M1_M2:.5f}, "
                    f"Max dM0/dt: {max_dM0_dt:.5f}, Max dM1/dt: {max_dM1_dt:.5f}, "
                    f"Max dM2/dt: {max_dM2_dt:.5f}, Max dD/dt: {max_dD_dt:.5f}"
                )




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