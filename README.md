# Percolation-diffusion model
MATLAB code to compute the diffusional re-equilibration of trace elements in a solid matrix of spherical mantle minerals percolated by a melt in a 1D column
    
The package (tested in Matlab R2021b and R2022b) includes:
- The main percolation-diffusion model (_mainPercolationDiffusion.m_)
- An input file (_input.xlsx_) which includes solid, liquid and normalization compositions, and a set of partition and diffusion coefficients
- A user file (_initializeAdhoc.m_) to load the inputs and set the different percolation-diffusion parameters 
- A _utils_ folder that contains all the required dependencies
- An _output_ folder where the models results are saved. This folder currently contains the pre-run model results necessary to reproduce the outputs shown in Tilhac et al. (currently under review)
- An additional code (_plotting.m_) to generate plots similar to the figures shown in Tilhac _et al_.

## Instructions & user changes

The model can be directly launched by running _mainPercolationDiffusion.m_. 
There are two main sections that contains parameters that can be safely changed by the user (marked as USER CHANGE):
- Firstly, the **Main parameters**, which includes P-T conditions, porosity, melt velocity, time and saving intervals, as well as the proportion of Eu<sup>2+</sup>/Eu<sup>3+</sup> (to investigate the diffusive fractionation of europium among the REE, as in Tilhac _et al._).
- Secondly, the **Mineral compositions**, which includes modal compositions, grain size and the activation of the P and T dependencies on diffusivities. The default mineral allocations is olivine (_Oli_), clinopyroxene (_Cpx_), orthopyroxene (_Opx_), garnet (_Grt_), spinel (_Spl_) and plagioclase (_Plg_), but only clinopyroxene is used in Tilhac _et al_. (see the **Reproducibility note** below).

The user can change the compositions and coefficients provided in the input file, but it is important that the partition and diffusion coefficients provided in the input file match the mineral compositions in the user file (_e.g._ if clinopyroxene, orthopyroxene and olivine are listed in the user file, their respective coefficients must be provided).

Running the model on a desktop computer using the current default settings takes about 1-10 minutes (depending, among others, on the number of time steps required). Do not force-change the number of time steps as needs to be calculated from the column length and melt velocity based on the number of nodes and particle spacing.

## Reproducibility note

Model outputs similar to the ones shown in Tilhac _et al_. can be reproduced by running _plotting.m_. To reproduce **Figure 1** (REE diagrams) and **Figure 2** (Eu anomalies _vs_ Eu contents), the code reads the output subfolder _pre-run_. Similar outputs can directly be obtained by running _mainPercolationDiffusion.m_ as provided, leaving _initializeAdhoc.m_ unchanged.

The _benchmark_Eu_ boolean allows to choose between two diffusivities for Eu<sup>2+</sup> based on the experimental diffusivities (Sneeringer _et al._, 1984) of Sr<sup>2+</sup> in either synthetic (1 - conservative choice used in Tilhac _et al._) or natural (0) diopside.

The _benchmark_2Cpx_ boolean allows to use two different grain-size populations for clinopyroxene (1 - used in Tilhac _et al._) or the default allocations (0).

## Credits

By Romain Tilhac, CSIC post-doctoral researcher at the  Instituto Andaluz de Ciencias de la Tierra (IACT), Granada, Spain.
Modified from the original MPMCRT code developed by Beñat Oliveira Bravo. <br />

Tilhac, R., Hidas, K., Oliveira, B., Garrido, C.J. Evidence of ghost plagioclase signature induced by kinetic fractionation of europium in the Earth’s mantle. Under review.

Sneeringer, M., Hart, S. R., Shimizu, N., 1984, Strontium and samarium diffusion in diopside: Geochimica et Cosmochimica Acta, v. 48, no. 8, p. 1589-1608.

Updated on December 9th, 2022 <br />
Contact: romain.tilhac@csic.es <br />
