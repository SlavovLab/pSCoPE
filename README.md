# **pSCoPE**

## Prioritized Single-Cell ProtEomics by Mass Spectrometry

<!--![GitHub release](https://img.shields.io/github/release/SlavovLab/DO-MS.svg)-->
![GitHub](https://img.shields.io/github/license/SlavovLab/DO-MS.svg)  |  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7498171.svg)](https://doi.org/10.5281/zenodo.7498171)

* [Website](https://scp.slavovlab.net/pSCoPE)
* [Running](#running)
* [Download data](https://scp.slavovlab.net/Huffman_et_al_2022)
* [Preprint](https://www.biorxiv.org/content/10.1101/2022.03.16.484655v1)



![workflow](Workflow_v7.png)

### Requirements

[MaxQuant.live (v. 2.1)](www.maxquant.live) contains the prioritization feature applied in the pSCoPE manuscript. This tool is compatible with Thermo orbitrap instruments running Tune v. 2.9 and Microsoft .Net v. 4.7.2. An automated system compatibility assessment is available on the linked MaxQuant.Live website.

The code blocks in this repository for figure generation and inclusion list generation have been tested on R >= 4.0.2, OSX 10.14 / Windows 10. R can be downloaded from the main [R Project page](https://www.r-project.org/) or downloaded with the [RStudio Application](https://www.rstudio.com/products/rstudio/download/).



### Running

#### Publication figure generation
The minimal data needed for figure generation is available by following the "Download data" link above and downloading the figure_data.zip directory.

All publication figures can be generated separately by running either the PDAC_BMDM_technical_figure_generation_script.rmd file or the BMDM_bio_figure_generation_script.rmd file, respectively.

#### Inclusion-list generation
A prioritized inclusion list can generated from DIA data searched by Spectronaut using the prioritizedInclusionListGen.rmd file. A sample set of Spectronaut search results to use can be found in the sampleFiles directory.

#### MaxQuant.Live log files
Information from the MaxQuant.Live log files can be extracted using the MaxQuantLive_Log_Stats_Script.rmd file.

#### Limited FASTA
The limited FASTA files used for searching prioritized single-cell experiments were generated using the limited_FASTA_generation_script.rmd file.

#### FASTA including MEROPS-annotated peptides
FASTA files containing MEROPS-annotated proteolytic cleavage products can be generated via the MEROPS_FASTA_generation_script.rmd file.



------------

The manuscript is freely available on bioRxiv: <!-- [Specht et al., 2019](https://www.biorxiv.org/content/10.1101/665307v3) -->

Contact the authors by email: [nslavov\{at\}northeastern.edu](mailto:nslavov@northeastern.edu).

### License

The pSCoPE code is distributed by an [MIT license](https://github.com/SlavovLab/DO-MS/blob/master/LICENSE).

### Contributing

Please feel free to contribute to this project by opening an issue or pull request.

<!--
### Data
All data used for the manuscript is available on [UCSD's MassIVE Repository](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=b15cafc7489147e99b93bd7c718388b2)
-->
