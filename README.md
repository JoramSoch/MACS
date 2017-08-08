# MACS

<h3>MACS – a new SPM toolbox for model assessment, comparison and selection</h3>

This toolbox unifies cross-validated Bayesian model selection (cvBMS) [1] and averaging (cvBMA) [2] as well as several other methods into a fully-functionining toolbox for Statistical Parametric Mapping (SPM) that supports SPM's batch mode to accelerate model assessment, comparison and selection (MACS) of general linear models (GLMs) applied to functional magnetic resonance imaging (fMRI) data.

This is <b>MACS V1.0</b>, also referred to as <b>MACS R2017a</b>, released on <b>24/05/2017</b>. The developers intend to immediately commit bug fixes to this repository and provide a general update two times a year. A toolbox manual and a preprint of the toolbox paper will soon be available. Example analyses and template batches can already be found in the subdirectories "MACS_Examples" and "MACS_Pipelines".

To install the toolbox, it has to be downloaded and placed as a subdirectory "MACS" into the SPM toolbox folder. Upon starting SPM, batch modules for toolbox features can be accessed by clicking "SPM -> Tools -> MACS Toolbox" in the SPM batch editor. MACS is optimized for SPM12, but also compatible with SPM8. With the adoption of MACS, cvBMS [3] and cvBMA [4] have become obsolete.

[1] http://www.sciencedirect.com/science/article/pii/S1053811916303615 <br>
[2] http://www.sciencedirect.com/science/article/pii/S105381191730527X <br>
[3] https://github.com/JoramSoch/cvBMS <br>
[4] https://github.com/JoramSoch/cvBMA <br>
