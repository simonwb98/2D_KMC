# src/kinetic_monte_carlo.py
import random 
import math

from sortedcontainers import SortedList

class RateContainer:
    def __init__(self):
        self.events = SortedList(key=lambda event: event[2])  # Sort by cumulative rate

    def add_event(self, monomer, action, rate):
        cumulative_rate = self.get_total_rate() + rate
        self.events.add((monomer, action, cumulative_rate))

    def get_total_rate(self):
        return self.events[-1][2] if self.events else 0

    def select_event(self, u):
        """
        Select an event probabilistically based on a random number u in [0, 1).
        """
        cumulative_rate = u * self.get_total_rate()
        for monomer, action, rate in self.events:
            if cumulative_rate <= rate:
                return monomer, action


def calculate_global_rate(lattice, monomers):
    """
    Calculate the total system rate and individual monomer rates.

    Args:
        lattice (Lattice): The lattice object.
        monomers (list): List of monomers.

    Returns:
        tuple: (R_total, events), where:
            - R_total is the total system rate.
            - events is a list of (monomer, action, rate) tuples.
    """
    R_total = 0
    events = []

    for monomer in monomers:
        # Calculate rates for all actions
        diffusion_rate = monomer.calculate_diffusion_rate(lattice)
        rotation_rate = monomer.calculate_rotation_rate(lattice)
        coupling_rate = monomer.calculate_coupling_rate(lattice)

        # Add to total rate and record events
        if diffusion_rate > 0:
            R_total += diffusion_rate
            events.append((monomer, "diffuse", diffusion_rate))
        if rotation_rate > 0:
            R_total += rotation_rate
            events.append((monomer, "rotate", rotation_rate))
        if coupling_rate > 0:
            R_total += coupling_rate
            events.append((monomer, "couple", coupling_rate))

    return R_total, events


def select_event(events, R_total):
    """
    Select an event probabilistically based on rates.

    Args:
        events (list): List of (monomer, action, rate) tuples.
        R_total (float): Total system rate.

    Returns:
        tuple: Selected (monomer, action).
    """
    u = random.random() * R_total
    cumulative_rate = 0

    for monomer, action, rate in events:
        cumulative_rate += rate
        if u <= cumulative_rate:
            return monomer, action

    raise RuntimeError("Event selection failed. Check rates and R_total.")

def advance_time(R_total):
    """
    Advance the system time.

    Args:
        R_total (float): Total system rate.

    Returns:
        float: Time step (Delta t).
    """
    u = random.random()
    return -math.log(u) / R_total

def perform_event(monomer, action, lattice):
    """
    Perform the selected event.

    Args:
        monomer (Monomer): The monomer involved in the action.
        action (str): The action to perform ("diffuse", "rotate", "couple").
        lattice (Lattice): The lattice object.
    """
    if action == "diffuse":
        monomer.diffuse(lattice)
    elif action == "rotate":
        monomer.rotate(lattice)
    elif action == "couple":
        monomer.couple(lattice)

def kmc_simulation(lattice, monomers, max_steps=1e5):
    """
    Perform the kinetic Monte Carlo simulation.

    Args:
        lattice (Lattice): The lattice object.
        monomers (list): List of monomers.
        max_steps (int): Maximum number of steps.

    Returns:
        float: Total simulation time.
    """
    time = 0

    for step in range(int(max_steps)):
        # Step 1: Calculate global rate and events
        R_total, events = calculate_global_rate(lattice, monomers)
        if R_total == 0:
            break  # No more possible events

        # Step 2: Select an event
        monomer, action = select_event(events, R_total)

        # Step 3: Advance time
        delta_t = advance_time(R_total)
        time += delta_t

        # Step 4: Perform the selected event
        perform_event(monomer, action, lattice)

    print(f"Simulation completed in {time:.2f} units.")
    return time

def kmc_simulation_optimized(lattice, monomers, max_steps=1e5):
    """
    Perform the optimized kMC simulation.
    """
    time = 0
    rate_container = RateContainer()

    # Precompute neighbors
    lattice.precompute_neighbors()

    # Initialize rates and events
    for monomer in monomers:
        monomer.update_rates(lattice)
        rate_container.add_event(monomer, "diffuse", monomer.cached_diffusion_rate)
        rate_container.add_event(monomer, "rotate", monomer.cached_rotation_rate)
        rate_container.add_event(monomer, "couple", monomer.cached_coupling_rate)

    for step in range(int(max_steps)):
        R_total = rate_container.get_total_rate()
        if R_total == 0:
            break  # No more events possible

        # Select an event
        u1 = random.random()
        monomer, action = rate_container.select_event(u1)

        # Perform the selected event
        perform_event(monomer, action, lattice)

        # Update affected rates
        monomer.update_rates(lattice)
        # Update the rate container dynamically...

        # Advance time
        u2 = random.random()
        delta_t = -math.log(u2) / R_total
        time += delta_t

    print(f"Optimized simulation completed in {time:.2f} units.")
    return time
