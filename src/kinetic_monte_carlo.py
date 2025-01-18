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
        print("Monomer diffused")
    elif action == "rotate":
        monomer.rotate(lattice)
        print("Monomer rotated")
    elif action == "couple":
        monomer.couple(lattice)
        print("Monomer coupled")

import random
import math

# Kinetic Monte Carlo simulation function
def kmc_simulation(lattice, monomers, max_steps=1e5):
    """
    Perform the kinetic Monte Carlo simulation.

    Args:
        lattice (Lattice): The lattice object.
        monomers (list): List of monomers.
        max_steps (int): Maximum number of simulation steps.

    Returns:
        float: Total simulation time.
    """
    time = 0

    for step in range(int(max_steps)):
        # Step 1: Calculate global rate and events
        R_total, events = calculate_global_rate(lattice, monomers)

        if R_total == 0:
            print("No more events possible. Simulation ended early.")
            break  # Exit if no events can occur

        # Step 2: Select an event probabilistically
        u1 = random.random()  # Uniform random number for event selection
        monomer, action = select_event(events, R_total, u1)

        # Step 3: Advance the simulation time
        u2 = random.random()  # Uniform random number for time advancement
        delta_t = -math.log(u2) / R_total
        time += delta_t

        # Step 4: Perform the selected event
        perform_event(monomer, action, lattice)

    print(f"Simulation completed in {time:.6f} time units after {step+1} steps.")
    return time

# Helper functions

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
        if not monomer.coupled:  # Skip coupled monomers
            diffusion_rate = monomer.calculate_diffusion_rate(lattice)
            rotation_rate = monomer.calculate_rotation_rate(lattice)
            coupling_rate = monomer.calculate_coupling_rate(lattice)

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


def select_event(events, R_total, u):
    """
    Select an event probabilistically based on rates.

    Args:
        events (list): List of (monomer, action, rate) tuples.
        R_total (float): Total system rate.
        u (float): Random number in [0, 1).

    Returns:
        tuple: Selected (monomer, action).
    """
    threshold = u * R_total
    cumulative_rate = 0

    for monomer, action, rate in events:
        cumulative_rate += rate
        if cumulative_rate >= threshold:
            return monomer, action

    raise RuntimeError("Event selection failed. Check rates and R_total.")


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

