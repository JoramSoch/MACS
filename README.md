# MACS

<a href="https://doi.org/10.5281/zenodo.2529423"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.2529423.svg" alt="DOI"></a>

<h3>MACS – a new SPM toolbox for model assessment, comparison and selection</h3>

This toolbox (pronounced as "Max") evaluates general linear models (GLMs) for functional magnetic resonance imaging (fMRI) data estimated in Statistical Parametric Mapping (SPM). MACS includes classical, information-theoretic and Bayesian methods of model assessment previously applied to GLMs for fMRI as well as recent methodological developments of model selection [1] and model averaging [2] in fMRI data analysis [3].

This is <b>MACS V1.3</b>, also referred to as <b>MACS R2018b</b>, released on <b>31/12/2018</b>. The developers intend to immediately commit bug fixes to this repository and provide a general update two times a year. A toolbox paper has been published in a peer-reviewed journal [3] and a toolbox manual is included in the repository [4].

To install the toolbox, it has to be downloaded and placed as a subdirectory "MACS" into the SPM toolbox folder. Upon starting SPM, batch modules for toolbox features can be accessed by clicking "SPM -> Tools -> MACS Toolbox" in the SPM batch editor [3, Fig. 3; 4, Fig. 1]. MACS is optimized for SPM12, but also compatible with SPM8.

The repository includes a number of sub-directories:
- `MACS_Examples`: SPM batch editor job files for example analyses from the toolbox paper [3, Sec. 4]
- `MACS_Pipelines`: SPM template batches/script for cvBMS [1], cvBMA [2] and model space definition
- `MACS_Extensions`: MATLAB scripts for toolbox extensions as described in the manual [4, Sec. 15]
- `MACS_Manual`: TEX and PDF file belonging to the latest version of the toolbox manual

[1] https://www.sciencedirect.com/science/article/pii/S1053811916303615 <br>
[2] https://www.sciencedirect.com/science/article/pii/S105381191730527X <br>
[3] https://www.sciencedirect.com/science/article/pii/S0165027018301468 <br>
[4] https://github.com/JoramSoch/MACS/blob/master/MACS_Manual/Manual.pdf <br>
