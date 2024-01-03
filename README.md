# Heterogeneous_diffusion

This is a MATLAB/Netlogo code for the paper "Beyond microtubules: The cellular environment at the endoplasmic reticulum attracts proteins to the nucleus, enabling nuclear transport" by Chae et al. (2023) submitted in Science.
The code simulates the protein diffusion in a heterogeneous intracellular environment with physical interaction with the ER. 

## Code Description
### Figure 2) Protein diffusion under the presence of obstacles.
1. Simulation_feedback.m
> Code for protein diffusion in an intracellular environment with 

### Figure 3) Effective PDE describing diffusion in the heterogeneous intracellular environment.
1. 
> PDE simulation code for Chapman's law.

2. diffusion2d_high_hetero.m
> MATLAB function code required for implementing pde_simulation_by_time.m
> This function calculates the flux and spatial derivative of PDE based on Chapman's law in the highly heterogeneous intracellular environment.


### Figure 4) ABM equivalent to effective PDE (Chapman's law)
1. hetero_diffusion_ellipse_v1.nlogo
> Simulates effective diffusion following Chapman's law in an elliptic-shaped cell.
2. hetero_diffusion_various_ER_v1.nlogo
> Simulates effective diffusion following Chapmna's law with a complex structured ER.

