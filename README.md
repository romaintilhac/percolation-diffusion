# Percolation-diffusion model
MATLAB code to compute the diffusional re-equilibration of trace elements in a solid matrix of spherical mantle minerals percolated by a melt in a 1D column
    
The package (tested in Matlab R2021b and R2022b) includes:
- The main percolation-diffusion model (mainPercolationDiffusion.m)
- An input file (input.xlsx) which includes solid, liquid and normalization compositions, and a set of partition and diffusion coefficients
- A user file (initializeAdhoc.m) to load the inputs and set the different percolation-diffusion parameters 
- A utils folder that contains all the required dependencies
- An output folder where the models results are saved. This folder currently contains the pre-run model results necessary to reproduce the outputs shown in Tilhac et al. (currently under review)
- An additional code (plotting.m) to generate plots similar to the figures shown in Tilhac et al.

## Instructions & user changes

The model can be directly launched by running mainPercolationDiffusion.m. 
There are two main sections that contains parameters that can be safely changed by the user (marked as USER CHANGE):
- Firstly, the Main parameters, which includes P-T conditions, porosity, melt velocity, time and saving intervals, as well as the proportion of Eu2+/Eu3+ (to investigate the diffusive fractionation of Eu as in Tilhac et al.).
- Secondly, the Mineral compositions, which includes modal compositions, grain size and the activation of the P and T dependencies on diffusivities.

The user can also change the compositions and coefficients provided in the input file. It is important that the partition and diffusion coefficients provided input file match the mineral compositions in the user file (e.g. if clinopyroxene, orthopyroxene and olivine are listed in the user file, they will need their respective coefficients to be provided). Note: the 'benchmark_2Cpx' boolean allows to use two different grain-size populations (using the two first allocations) for cpx, as in Tilhac et al. (see the Reproducibility note below).

Running the model on a desktop computer using the current default settings takes about 1-10 minutes (depending, among others, on the number of time steps required). Do not force-change the number of time steps as needs to be calculated from the column length and melt velocity based on the number of nodes and particle spacing.

## Reproducibility note

Model outputs similar to the ones shown in Tilhac et al. can be reproduced by running plotting.m.
- To reproduce Figures 1 (REE diagrams) and 2 (Eu anomalies vs Eu content), the code reads the output subfolder 'Oran_2cpx'. Similar outputs can directly be obtained by running 'mainPercolationDiffusion.m' as provided, leaving initializeAdhoc.m unchanged.

Note: the 'benchmark_Eu' boolean allows to choose between two diffusivities for Eu2+ based on the experimental diffusivities (Sneeringer et al., 1984) of Sr2+ in either synthetic (1 - conservative choice used in Tilhac et al.) or natural (0) diopside.

## Credits

Modified from the original MPMCRT code developed by Beñat Oliveira Bravo <br />
By Romain Tilhac, CSIC post-doctoral researcher at the  Instituto Andaluz de Ciencias de la Tierra (IACT), Granada, Spain.

Tilhac, R., Hidas, K., Oliveira, B., Garrido, C.J. Evidence of ghost plagioclase signature induced by kinetic fractionation of europium in the Earth’s mantle. Under review.

Sneeringer, M., Hart, S. R., Shimizu, N., 1984, Strontium and samarium diffusion in diopside: Geochimica et Cosmochimica Acta, v. 48, no. 8, p. 1589-1608.

Updated on December 9th, 2022 <br />
Contact: romain.tilhac@csic.es <br />
