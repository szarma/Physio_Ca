# Physio_Ca


A toolbox to analyze and interact with Ca imaging data, developed within the Cell and Tissue Networks research group led by [prof. Marjan Slak-Rupnik](https://www.meduniwien.ac.at/web/index.php?id=688&res_id=37&name=Marjan_Slak%20Rupnik) at the Medical University of Vienna. 

https://user-images.githubusercontent.com/2512087/162633046-b26d7c49-3501-4e78-9a33-433119157537.mp4

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

## Documentation

For brief instructions on how to process a recording, have a look at [Step 0 of processing](notebooks/Tutorials/Processing_Tutorial.html). 

This will get your recording motion corrected and segmented, potentially with hundreds of ROIs. You will then want to clean them up a bit ([Step 1](notebooks/Tutorials/Step1_roi_massages.html)). 

Then, you may want to detect all transients (i.e. events) in the traces ([Step 2](notebooks/Tutorials/Step2_rois2events.html)), and do some basic exploratory analysis and visualizations of the events ([Step 3](notebooks/Tutorials/Step3_event_visualization_and_analysis.html)).

For a more didactical introduction into filtering and event detection, you may want to check [this page](docs/events.html).
And for a more indepth dive into more technical details, have a look [here](docs/matmet.pdf).

For answering actual scientific questions, we need to pool different experiments. In a simple case, where each experiment is to be treated as independent replicate, you can follow [this example analysis](notebooks/Tutorials/Step4_experiment_pooling.html).

If the experiments differ in some important variable (e.g. some mice are treated, and are not), and we are interested in quantifying the difference between them (i.e. effect of the variable of interest), you can follow [this example analysis](notebooks/Tutorials/Step4_pooling_multi_legs_and_voi.html).

## Features

One of our most used tools is the _''roi examiner''_ for live interaction with the ROIs and their traces within a jupyter notebook.
![examiner_demo](https://user-images.githubusercontent.com/2512087/162623035-c054b171-c222-47b0-905e-6f91fcb0caab.gif)

We have also build a dashboard for fast intreaction with our storage filesystem. Given a folder, it finds all processed recordings in it and its subfolders, collects metadata and presents it in a table form. It further enables entering of the experimental protocol, and additional data, which are then also searchable. It also provides a link to an automaticaly generated notebook for a brief glimpse into actual results of an experiment.





## Docker
We are working on a docker container with all important bits of our toolbox installed. Check the [docker directory](docker/) for more details.
