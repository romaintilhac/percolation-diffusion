# Percolation-diffusion model [![DOI](https://zenodo.org/badge/573516734.svg)](https://zenodo.org/badge/latestdoi/573516734)
MATLAB code to compute the diffusional re-equilibration of trace elements (_e.g._ REE) in a solid matrix of spherical mantle minerals percolated by a melt in a 1D column
    
The package (tested in Matlab R2021b and R2022b) includes:
- The main percolation-diffusion model (_mainPercolationDiffusion.m_)
- An input file (_input.xlsx_) which includes solid, liquid and normalization compositions, and a set of partition and diffusion coefficients for the main mantle minerals
- A user file (_initializeAdhoc.m_) to load the inputs and set the different percolation-diffusion parameters 
- A _utils_ folder that contains all the required dependencies
- An _output_ folder where the models results are saved. This folder currently contains the pre-run model results necessary to reproduce the outputs shown in [Tilhac _et al._ (2023)](https://www.nature.com/articles/s41467-023-36753-0).
- An additional code (_plotting.m_) to generate plots similar to the figures shown in [Tilhac _et al._ (2023)](https://www.nature.com/articles/s41467-023-36753-0).

## Instructions & user changes

The model can be directly launched by running _mainPercolationDiffusion.m_ (only the size _y2_ [m] of the column may be changed in this file).
There are two main sections that contains parameters that can be safely changed by the user in _initializeAdhoc.m_ (marked as USER CHANGE):
- Firstly, the **Main parameters**, which include P-T conditions, porosity, melt velocity, time and saving intervals, as well as the proportion of Eu<sup>2+</sup>/Eu<sup>3+</sup> (to investigate the diffusive fractionation of europium among the REE, as in [Tilhac _et al._, 2023](https://www.nature.com/articles/s41467-023-36753-0)).
- Secondly, the **Mineral parameters**, which include modal proportions, grain size and the activation of the P and T dependencies on diffusivities. The default mineral allocations are olivine (_Oli_), clinopyroxene (_Cpx_), orthopyroxene (_Opx_), garnet (_Grt_), spinel (_Spl_) and plagioclase (_Plg_).

The partition and diffusion coefficients provided in the input file are taken from [Oliveira _et al._ (2020)](https://doi.org/10.1093/petrology/egaa067). They can be changed by the user but it is important that they match the mineral parameters in the user file (_e.g._ if clinopyroxene, orthopyroxene and olivine are listed in the user file, their respective coefficients must be provided). The input compositions given as example correspond to the data used in [Tilhac _et al._ (2023)](https://www.nature.com/articles/s41467-023-36753-0) (see **Reproducibility note** below).

Running the model on a desktop computer using the current default settings takes about 1-10 minutes (depending, among others, on the number of time steps required). Do not manually change the number of time steps as it needs to be calculated from the column length and melt velocity based on the number of nodes and particle spacing.

## Reproducibility note

Model outputs as shown in [Tilhac _et al._ (2023)](https://www.nature.com/articles/s41467-023-36753-0) can be reproduced by directly running _plotting.m_. To reproduce **Figure 1** (REE patterns) and **Figure 2** (Eu anomalies _vs_ Eu contents), the code reads the output subfolder _pre-run_. Similar outputs can also directly be obtained by running _mainPercolationDiffusion.m_ as provided, leaving _initializeAdhoc.m_ unchanged.

The _benchmark_Eu_ boolean allows to choose between two diffusivities for Eu<sup>2+</sup> based on the experimental diffusivities of Sr<sup>2+</sup> in either synthetic (1 - conservative choice used in [Tilhac _et al._, 2023](https://www.nature.com/articles/s41467-023-36753-0)) or natural (0) diopside.

The _benchmark_2Cpx_ boolean allows to use two different grain-size populations for clinopyroxene only (1 - used in [Tilhac _et al._, 2023](https://www.nature.com/articles/s41467-023-36753-0)) or the default mineral allocations (0).

## Credits

By Romain Tilhac, CSIC post-doctoral researcher at the  Instituto Andaluz de Ciencias de la Tierra (IACT), Granada, Spain. <br />
Modified from the original MPMCRT code developed by Beñat Oliveira Bravo ([Oliveira _et al.,_ 2020](https://doi.org/10.1093/petrology/egaa067)). <br /> 
Updated on December 11th, 2022 <br />
Contact: romain.tilhac@csic.es

### References

Tilhac, R., Hidas, K., Oliveira, B., Garrido, C.J. 2023. Evidence of ghost plagioclase signature induced by kinetic fractionation of europium in the Earth’s mantle, Nature Communications, Volume 14, 1099, https://www.nature.com/articles/s41467-023-36753-0

Oliveira, B., Afonso, J.C., Tilhac, R. 2020. A Disequilibrium Reactive Transport Model for Mantle Magmatism, Journal of Petrology, Volume 61, Issue 9, egaa067, https://doi.org/10.1093/petrology/egaa067

