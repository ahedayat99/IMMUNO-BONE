"""
ElementAgent: a mixed-composition population parcel (Methods 5.1, 5.5).

An agent is not a single cell. Each parcel resides in one finite element and carries continuous
amounts of five cell populations together with the local cytokine state it senses:

    state = [PMN, M0, M1, M2, c1, c2, c3, MSC, c4]
             cells (cell-equivalents)  ^ cytokines (copy of the element fields: TNF-a, IL-10, TGF-b, VEGF-A)

The cell amounts are advanced by the element reaction step (DomainModel.react, scripts/reaction.py),
where macrophage polarization moves amount between the M0, M1 and M2 compartments at rates given by
the Hill functions of Trejo et al. 2019 (the same functional forms are used in COMMBINI, Borgiani et al. 2023).
ElementAgent.step() applies the phenotype-budding rules; migration is DomainModel.migrate_parcels().
"""

import numpy as np
from mesa import Agent

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
        # Totals are computed once per step by the model and updated as buds are created,
        # instead of re-scanning every parcel on every bud (same cap, linear cost).
        totals = self.model.population_totals
        if totals[state_index] + amount_to_add > self.MAX_POPULATION:
            return False
        totals[state_index] += amount_to_add
        return True

    def step(self):
        """
        Parcel events after the element reaction and transport steps (Methods 5.5, Box 1 step 4a):
        phenotype budding. The parcel's cell amounts were already advanced by DomainModel.react();
        its cytokine slots (4, 5, 6, 8) hold a copy of the element fields.
        """
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
                    self.model.remove_parcel(self)
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
                    self.model.remove_parcel(self)
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
                    self.model.remove_parcel(self)
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
