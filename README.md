# Pricing and Energy Trading in Peer-to-peer Zero Marginal-cost Microgrids

This repository provides the models used in the paper: J. T. Lee, R. Henriquez-Auba, B. K. Poolla and D. S. Callaway, "Pricing and Energy Trading in Peer-to-peer Zero Marginal-cost Microgrids", submitted to IEEE Transactions on Smart Grid.

Please note that this repository already includes a copy of the repository [Computational Experiments in MATLAB](https://github.com/leejt489/computational-experiment-matlab) to simplify the simulation of the cases presented in the paper.

## Run Simulation Cases

To run the simulations in MATLAB presented in the paper in Section IV, the user should first clone the repository to its local computer. Afterwards, run the file `init.m` to properly set-up the path in MATLAB.

### A) Effect of parameters γ and δ₀ on convergence.

To run this simulation, open the file `analyzeParameterExperiment.m` and run it. Please note that this will run 100 trials for each of the 380 pairs of (γ, δ₀) and the simulation will take around 8 hours.

### B) Performance for unproven cases

To be added soon.

## Algorithm implementations

Algorithms' implementation are provided the `code` folder with comments. `Bid_process.m` implements Algorithm 1 of the paper. `fit_utility.m` is used to construct the utility functions based on the price, quantity and elasiticity. `opt_centralized.m` provides the implementation of the centralized formulation presented in Section II. `opt_pi.m` and `opt_q.m` solves the pi-agent and q-agent optimization problems, respectively, required in Algorithm 1.