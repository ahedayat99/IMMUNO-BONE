"""
ElementAgent class for element-based agent simulation - OPTIMIZED VERSION
"""

import numpy as np
from mesa import Agent
from scipy.integrate import solve_ivp
from bone_healing_model_optimized import full_bone_healing_model

# Import at module level instead of in step() method
from endothelial_cell_agent import EndothelialCellAgent


# Define smooth_transition at module level to avoid redefining
def smooth_transition(current_value, target_value, rate=0.1):
    """Helper function for smooth state transitions."""
    return current_value + rate * (target_value - current_value)


class ElementAgent(Agent):
    # Population caps for each cell type
    MAX_POPULATION = 10000

    def __init__(self, unique_id, model, initial_conditions, params, element_id, centroid):
        """
        Initialize an agent associated with an FE element.
        :param unique_id: Unique identifier for the agent.
        :param model: The model to which the agent belongs.
        :param initial_conditions: Initial state variables for the agent.
        :param params: Model parameters for the bone healing equations.
        :param element_id: The ID of the element this agent belongs to.
        :param centroid: The centroid of the element in the mesh.
        """
        super().__init__(unique_id, model)
        self.state = np.array(initial_conditions)
        self.params = params
        self.element_id = element_id
        self.centroid = centroid

        # Cache agent type to avoid repeated string operations
        self.agent_type = unique_id.split('-')[0] if '-' in unique_id else unique_id

    def _get_total_population(self, state_index):
        """
        Get total population of a specific cell type across all agents.
        :param state_index: Index in state array (0=PMN, 1=M0, 2=M1, 3=M2, 7=MSC)
        :return: Total count of that cell type
        """
        total = 0
        for agent in self.model.schedule.agents:
            if isinstance(agent, ElementAgent):
                total += agent.state[state_index]
        return total

    def _can_spawn(self, state_index, amount_to_add):
        """
        Check if spawning would exceed population cap.
        :param state_index: Index in state array
        :param amount_to_add: Amount that would be added
        :return: True if spawning is allowed, False otherwise
        """
        current_total = self._get_total_population(state_index)
        return (current_total + amount_to_add) <= self.MAX_POPULATION

    def step(self):
        """
        Perform one simulation step for this agent, updating its state.
        """
        t_span = (0, 0.1)

        # Use the element ID to get the debris concentration
        D = self.model.debris_field[self.element_id]

        # Get EC count from pre-calculated cache (if available) or calculate
        if hasattr(self.model, 'element_EC_counts'):
            EC = self.model.element_EC_counts.get(self.element_id, 0)
        else:
            # Fallback to old method if cache not available
            EC = sum(1 for a in self.model.element_agents[self.element_id] if isinstance(a, EndothelialCellAgent))

        PO2 = self.model.oxygen_field[self.element_id]

        # Solve the ODE system
        solution = solve_ivp(
            lambda t, y: full_bone_healing_model(t, y, self.params, D, EC, PO2),
            t_span,
            self.state,
            method="RK45",
            atol=1e-6,
            rtol=1e-4,
            vectorized=True
        )

        # Update the agent's state with the result of the ODE integration
        updated_state = solution.y[:, -1]

        # Use cached agent type instead of string comparison
        if self.agent_type == "MSC":
            # Only update PMN, c1, c2, c3, c4, and Cm
            self.state[0] = updated_state[0]  # PMN
            self.state[4] = updated_state[4]  # c1
            self.state[5] = updated_state[5]  # c2
            self.state[6] = updated_state[6]  # c3
            self.state[8] = updated_state[8]  # c4
            self.state[7] = updated_state[7]  # Cm
        else:
            # For non-MSC agents, update all state variables
            self.state = updated_state

        PMN, M0, M1, M2, c1, c2, c3, Cm, c4 = self.state

        # Use cached agent type for all transitions
        # PMN to M0 transition
        if self.agent_type == "PMN":
            if self.state[1] > 1:
                transition_rate = 0.02
                M0_increase = smooth_transition(0, min(self.state[1], self.state[0]), rate=transition_rate)
                M0_spawn_amount = M0_increase * 10

                # Check if spawning would exceed cap
                if self._can_spawn(1, M0_spawn_amount):  # state_index 1 = M0
                    self.state[1] -= M0_increase
                    self.state[0] -= M0_increase * 10
                    self.model.spawn_agent(self.element_id, [0, M0_spawn_amount, 0, 0, c1, c2, c3, Cm, c4])

            if self.state[2] > 1:
                transition_rate = 0.02
                M0_increase = smooth_transition(0, min(self.state[2], self.state[0]), rate=transition_rate)
                M1_spawn_amount = M0_increase * 10

                # Check if spawning would exceed cap
                if self._can_spawn(2, M1_spawn_amount):  # state_index 2 = M1
                    self.state[2] -= M0_increase
                    self.state[0] -= M0_increase * 10
                    self.model.spawn_agent(self.element_id, [0, 0, M1_spawn_amount, 0, c1, c2, c3, Cm, c4])

        # M0 to M1 transition
        elif self.agent_type == "M0":
            if self.state[2] > 1:
                transition_rate = 0.03
                M1_increase = smooth_transition(0, min(self.state[2], self.state[1]), rate=transition_rate)
                M1_spawn_amount = M1_increase * 15
                M2_spawn_amount = M1_increase * 6

                # Check if spawning would exceed caps for M1 or M2
                if self._can_spawn(2, M1_spawn_amount) and self._can_spawn(3, M2_spawn_amount):
                    self.state[2] -= M1_increase
                    self.state[1] -= M1_increase * 100
                    self.model.spawn_agent(self.element_id, [0, 0, M1_spawn_amount, M2_spawn_amount, c1, c2, c3, M1_increase, c4])
                    if self in self.model.elements[self.element_id]:
                        self.model.elements[self.element_id].remove(self)
                    self.model.schedule.remove(self)
                    return

            # M0 to M2 transition (if direct transition allowed)
            if self.state[3] > 1:
                transition_rate = 0.05
                M2_increase = smooth_transition(0, min(self.state[3], self.state[1]), rate=transition_rate)
                M2_spawn_amount = M2_increase * 10

                # Check if spawning would exceed cap for M2
                if self._can_spawn(3, M2_spawn_amount):
                    self.state[3] -= M2_increase
                    self.state[1] -= M2_increase * 10
                    self.model.spawn_agent(self.element_id, [0, 0, 0, M2_spawn_amount, c1, c2, c3, M2_increase, c4])
                    if self in self.model.elements[self.element_id]:
                        self.model.elements[self.element_id].remove(self)
                    self.model.schedule.remove(self)
                    return

        # M1 to M2 transition
        elif self.agent_type == "M1":
            if self.state[3] > 1:
                transition_rate = 0.05
                M2_increase = smooth_transition(0, min(self.state[3], self.state[2]), rate=transition_rate)
                M2_spawn_amount = M2_increase * 10

                # Check if spawning would exceed cap for M2
                if self._can_spawn(3, M2_spawn_amount):
                    self.state[3] -= M2_increase
                    self.state[2] -= M2_increase * 2
                    self.model.spawn_agent(self.element_id, [0, 0, 0, M2_spawn_amount, c1, c2, c3, M2_increase, c4])
                    if self in self.model.elements[self.element_id]:
                        self.model.elements[self.element_id].remove(self)
                    self.model.schedule.remove(self)
                    return

        # Spawning logic for MSC agents with smooth transition
        elif self.agent_type == "MSC":
            if self.state[7] > 6:
                transition_rate = 0.05
                Cm_increase = smooth_transition(0, self.state[7] - 6, rate=transition_rate)
                MSC_spawn_amount = 10 * Cm_increase

                # Check if spawning would exceed cap for MSC
                if self._can_spawn(7, MSC_spawn_amount):  # state_index 7 = MSC (Cm)
                    self.state[7] -= 30*Cm_increase
                    self.model.spawn_agent(self.element_id, [0, 0, 0, 0, self.state[4], self.state[5], self.state[6], MSC_spawn_amount, self.state[8]])

        # Note: The debug calculations at the end of the original step() were removed
        # as they duplicate work done in domain_model.step()

    def migrate(self, migration_coefficient, random_migration_probability=0.1):
        """
        Move the agent based on the gradient of cytokines and debris between elements.
        MSC agents migrate toward high concentrations of c3 specifically.
        Include a small random chance to move to a neighboring element after migration.
        :param migration_coefficient: The proportionality factor for migration based on the gradient.
        :param random_migration_probability: Probability of performing a random migration.
        """
        # Current element and concentration
        current_element = self.element_id

        if self.agent_type == "MSC":
            # MSC agents migrate toward high concentrations of c3 only
            neighbors = self.model.get_neighbors(current_element)
            best_gradient = float('-inf')
            best_neighbor = None

            for neighbor in neighbors:
                if len(self.model.element_agents[neighbor]) < 3:
                    gradient = self.model.cytokine_fields["c3"][neighbor] - self.model.cytokine_fields["c3"][current_element]
                    if gradient > best_gradient:
                        best_gradient = gradient
                        best_neighbor = neighbor

            # If a valid migration target is found, move to it
            if best_neighbor is not None and best_gradient > 0:
                self._move_to_element(best_neighbor)
                return
        else:
            # Original migration logic for non-MSC agents
            current_concentration = self.model.debris_field[current_element] + \
                                    sum(self.model.cytokine_fields[c][current_element] for c in ['c1', 'c2', 'c3', 'c4'])
            neighbors = self.model.get_neighbors(current_element)

            neighbor_gradients = []
            for neighbor in neighbors:
                if len(self.model.element_agents[neighbor]) < 3:
                    neighbor_concentration = self.model.debris_field[neighbor] + \
                                            sum(self.model.cytokine_fields[c][neighbor] for c in ['c1', 'c2', 'c3', 'c4'])
                    gradient = neighbor_concentration - current_concentration
                    neighbor_gradients.append((neighbor, gradient))

            neighbor_gradients.sort(key=lambda x: x[1], reverse=True)

            for neighbor, gradient in neighbor_gradients:
                if gradient > 0:
                    self._move_to_element(neighbor)
                    self._random_migration(random_migration_probability)
                    return

            # Expand to second-level neighbors if no suitable immediate neighbor is found
            visited = set(neighbors)
            for neighbor, _ in neighbor_gradients:
                secondary_neighbors = self.model.get_neighbors(neighbor)
                for secondary_neighbor in secondary_neighbors:
                    if secondary_neighbor not in visited and len(self.model.element_agents[secondary_neighbor]) < 3:
                        secondary_concentration = self.model.debris_field[secondary_neighbor] + \
                                                sum(self.model.cytokine_fields[c][secondary_neighbor] for c in ['c1', 'c2', 'c3', 'c4'])
                        gradient = secondary_concentration - current_concentration
                        if gradient > 0:
                            self._move_to_element(secondary_neighbor)
                            self._random_migration(random_migration_probability)
                            return
                        visited.add(secondary_neighbor)

        # No suitable element found
        self._random_migration(random_migration_probability)

    def _move_to_element(self, target_element):
        """
        Helper function to handle agent migration to a new element.
        Ensures the element_agents list is updated properly.
        """
        # Remove agent from the current element
        self.model.element_agents[self.element_id].remove(self)

        # Update agent's position
        self.element_id = target_element
        self.centroid = self.model.element_centroids[target_element]

        # Add agent to the new element
        self.model.element_agents[target_element].append(self)

    def _random_migration(self, probability):
        """
        Helper function to perform a random migration to a neighboring element with a given probability.
        :param probability: The chance of performing a random migration.
        """
        if self.random.random() < probability:
            neighbors = self.model.get_neighbors(self.element_id)
            valid_neighbors = [neighbor for neighbor in neighbors if len(self.model.element_agents[neighbor]) < 3]

            if valid_neighbors:
                random_choice = self.random.choice(valid_neighbors)
                self._move_to_element(random_choice)
