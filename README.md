# Optimal_DG_Placement_ZBLF_Simulation
Simulation of "Optimum DG Placement for Known Power Injection From Utility/Substation by a Novel Zero Bus Load Flow Approach"
This repository includes the MATLAB implementation of the simulation based on the paper that introduces a two-part methodology: ZBLF (Zero Bus Load Flow) and SOS (Symbiotic Organisms Search) for optimal DG (Distributed Generation) placement.

Part 1: ZBLF Method
Objective: Find the optimal location and capacity of a single DG to minimize active power losses under known upstream injection and varying load scenarios.

Networks Studied: IEEE 33-bus and 69-bus distribution systems.

Method: Uses a modified Newton-Raphson algorithm where the utility bus serves as the voltage angle reference and the DG bus acts as a slack bus. The Jacobian matrix is adjusted accordingly.

Scenario Simulated: Zero power injection from the utility with peak load profile.

Optimization: Full-space search for optimal bus placement.

Results: The optimal DG location, power output, and reduced losses are reported.

Part 2: SOS Algorithm
Objective: Determine the optimal placement and power sharing of two DGs to further reduce system losses.

Method: Builds upon the ZBLF structure. One DG is fixed at the optimal ZBLF location, and the second is placed using sensitivity analysis. Power sharing is optimized using the SOS algorithm.

Network: IEEE 33-bus system.

Scenario: Assumes equal power factor for load and DG.

All simulations were implemented in MATLAB using M-files,
