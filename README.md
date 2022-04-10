# Physio_Ca

A toolbox to analyze and interact with Ca imaging data.
It developed within the Cell and Tissue Networks research group led by [prof. Marjan Slak-Rupnik](https://www.meduniwien.ac.at/web/index.php?id=688&res_id=37&name=Marjan_Slak%20Rupnik) at the Medical University of Vienna. 

A typical experiment involving imaging of pancreatic slices in our lab concerns a single field of view
showing up to hundreds of cells, in a recording of at least several, often dozens, gigabytes.
Current tools (i.e. ImageJ) rely on loading the recording, or its part, into memory, for viewing, analysis, and processing.
It also requires laborious and long human engagement.
We have developed a set of interdependent tools to automatize as much as possible the analysis pipeline. 

The main elements of our pipeline are the following:
 - Automatic detection of regions of interest (ROIs);
 - Transformation of ROI time traces into standard score ("z-score") and correction for filtering distortion;
 - Representation of the phenotype for each trace (height, auc, or halfwidth statistics, event rate...).
 - Quantification of the effect magnitude for acute pharmacological manipulations, and/or different experimental branches (different mice, or diet, of genetic manipulations).

![pipeline](https://user-images.githubusercontent.com/2512087/162617713-efd571a5-784e-4b2c-99ee-663f25457527.png)



And [here](docker/deployment.md) is page dedicated to docker deployment.
