# Percolation-diffusion model
MATLAB code to compute the diffusional re-equilibration of REE in a solid matrix of spherical mantle minerals percolated by a melt in a 1D column
    
The package (tested in Matlab R2021b) includes:
- The main percolation-diffusion model ('mainPercolationDiffusion.m')
- An input file (input.xlsx) which includes solid, liquid and normalization compositions, and a set of partition and diffusion coefficients. 
- A user file ('initializeAdhoc.m') to load the inputs and initialize the different percolation-diffusion parameters 
- A utils folder that containts ll the required dependencies
- An output folder that containts pre-run results similar to the model results shown in Tilhac et al. (currently under review).
- An additional code ('plotting.m')to generate plots similar to the figures shown in Tilhac et al. (currently under review).

The model can be directly launched by running 'mainPercolationDiffusion.m'. The user can change the input file (as long as the user file is updated accordingly). Running the model on a desktop computer using the current default settings takes about 1-10 minutes (depending on the number of timesteps required by the percolation settings).

Figures 1 (REE diagrams), 2 (Eu anomalies vs Eu content) and 3b (bulk Eu anomalies vs time) from Tilhac et al. can be reproduced by running 'plotting.m' directly using the pre-run results in the output folder. 

Modified from the original MPMCRT code developped by Beñat Oliveira Bravo <br />
By Romain Tilhac, Decembre 2022 <br />
Contact: romain.tilhac@csic.es <br />
