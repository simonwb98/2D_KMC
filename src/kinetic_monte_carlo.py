# src/kinetic_monte_carlo.py
import random 
import math


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

