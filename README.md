# MACS

<a href="https://zenodo.org/badge/latestdoi/92269905"><img src="https://zenodo.org/badge/92269905.svg" alt="DOI"></a>

<h3>MACS – a new SPM toolbox for model assessment, comparison and selection</h3>

This toolbox unifies cross-validated Bayesian model selection (cvBMS) [1] and averaging (cvBMA) [2] as well as several other methods into a fully-functionining toolbox for Statistical Parametric Mapping (SPM) that supports SPM's batch mode to accelerate model assessment, comparison and selection (MACS) of general linear models (GLMs) applied to functional magnetic resonance imaging (fMRI) data.

This is <b>MACS V1.2</b>, also referred to as <b>MACS R2018a</b>, released on <b>30/06/2018</b>. The developers intend to immediately commit bug fixes to this repository and provide a general update two times a year. A toolbox paper has been published in a peer-reviewed journal [3] and a toolbox manual is included in the repository [4].

To install the toolbox, it has to be downloaded and placed as a subdirectory "MACS" into the SPM toolbox folder. Upon starting SPM, batch modules for toolbox features can be accessed by clicking "SPM -> Tools -> MACS Toolbox" in the SPM batch editor [3, Fig. 3; 4, Fig. 1]. MACS is optimized for SPM12, but also compatible with SPM8. With the adoption of MACS, toolkits for cvBMS and cvBMA have become obsolete.

The repository includes a number of sub-directories:
- `MACS_Examples`: SPM batch editor job files for example analyses from the toolbox paper [3, Sec. 4]
- `MACS_Pipelines`: SPM template batches/script for cvBMS, cvBMA and model space definition
- `MACS_Extensions`: MATLAB scripts for toolbox extensions as described in the manual [4, Sec. 15]
- `MACS_Manual`: TEX and PDF file belonging to the latest version of the toolbox manual

[1] https://www.sciencedirect.com/science/article/pii/S1053811916303615 <br>
[2] https://www.sciencedirect.com/science/article/pii/S105381191730527X <br>
[3] https://www.sciencedirect.com/science/article/pii/S0165027018301468 <br>
[4] https://github.com/JoramSoch/MACS/blob/master/MACS_Manual/Manual.pdf <br>
