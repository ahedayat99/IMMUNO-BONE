"""
EndothelialCellAgent class for endothelial cell simulation.
"""

import numpy as np
from mesa import Agent


class EndothelialCellAgent(Agent):
    def __init__(self, unique_id, model, element_id, centroid):
        super().__init__(unique_id, model)
        self.element_id = element_id
        self.centroid = centroid
        self.sprout_path = [centroid]
        self.persistence = np.random.rand()
        self.prev_direction = np.array([1.0, 0.0])

    def step(self):
        neighbors = self.model.get_neighbors(self.element_id)
        current_pos = np.array(self.centroid)

        alpha = 2.0
        beta = 3.0

        weights = []
        direction_vectors = []
        candidate_neighbors = []

        for neighbor in neighbors:
            neighbor_pos = np.array(self.model.element_centroids[neighbor])
            direction = neighbor_pos - current_pos
            distance = np.linalg.norm(direction)

            if distance == 0:
                continue

            unit_direction = direction / distance
            vegf_diff = self.model.cytokine_fields["c4"][neighbor] - self.model.cytokine_fields["c4"][self.element_id]
            angle_alignment = np.dot(unit_direction, self.prev_direction)

            score = alpha * vegf_diff + beta * angle_alignment
            weights.append(np.exp(score))
            direction_vectors.append(unit_direction)
            candidate_neighbors.append(neighbor)

        total_weight = sum(weights)
        if total_weight == 0:
            return

        probabilities = [w / total_weight for w in weights]
        selected_idx = np.random.choice(len(candidate_neighbors), p=probabilities)

        chosen_element = candidate_neighbors[selected_idx]
        chosen_vector = direction_vectors[selected_idx]

        # Move and update trail
        self.model.move_agent_to(self, chosen_element)
        self.prev_direction = chosen_vector

        # Now tracking maturity and VEGF history too
        new_segment = {
            "path": self.sprout_path.copy(),
            "maturity": "immature",
            "age": 0,
            "last_high_VEGF_time": self.model.time
        }
        self.model.vessel_segments.append(new_segment)
        self.sprout_path.append(self.centroid)

        if self.element_id not in self.model.vessel_element_ids:
            self.model.vessel_element_ids.add(self.element_id)

        # 1. Stochastic branching
        if np.random.rand() < 0.05 and len(self.sprout_path) > 2:
            branch_target = self.model.find_branch_location(self.element_id)
            if branch_target:
                new_id = f"{self.unique_id}_branch_{self.model.time}"
                new_agent = EndothelialCellAgent(
                    new_id,
                    self.model,
                    branch_target,
                    self.model.element_centroids[branch_target]
                )
                self.model.schedule.add(new_agent)
                self.model.element_agents[branch_target].append(new_agent)
                print(f"{self.unique_id} branched at {self.element_id} → {new_id}")

        # 2. Anastomosis detection
        for other in self.model.schedule.agents:
            if isinstance(other, EndothelialCellAgent) and other != self:
                dist = np.linalg.norm(np.array(self.centroid) - np.array(other.centroid))
                if dist < 0.3:
                    print(f"Anastomosis: {self.unique_id} met {other.unique_id} at {self.centroid}")
                    self.model.link_sprouts(self, other)
                    break
