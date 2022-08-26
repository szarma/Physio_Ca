=========
Physio_Ca
=========

|PyPI Version badge| |PyPI Python Version badge| |Docker Image Version badge| |Docker Image Size badge|

A toolbox to analyze and interact with Ca imaging data, developed within the Cell and Tissue Networks research group led by `prof. Marjan Slak-Rupnik <https://www.meduniwien.ac.at/web/index.php?id=688&res_id=37&name=Marjan_Slak%20Rupnik>`_ at the Medical University of Vienna. 

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

.. image:: https://user-images.githubusercontent.com/2512087/162617713-efd571a5-784e-4b2c-99ee-663f25457527.png


Documentation
=============

The usage of the framework in practical terms is documented in the original repository at `https://github.com/szarma/Physio_Ca/ <https://github.com/szarma/Physio_Ca/>`_


Features
========

One of our most used tools is the _''roi examiner''_ for live interaction with the ROIs and their traces within a jupyter notebook.

.. image:: https://user-images.githubusercontent.com/2512087/162623035-c054b171-c222-47b0-905e-6f91fcb0caab.gif

We have a similar app to examine line scans.

.. image:: https://user-images.githubusercontent.com/2512087/162633612-ad71e643-14bb-4e62-b0f0-21188ec4c10c.gif

For examine detected events, one at a time, we also have a app.

.. image:: https://user-images.githubusercontent.com/2512087/162635307-6dea02ec-c56f-41ed-a275-efee595c1b9a.gif

We have also built a dashboard for fast intreaction with our storage filesystem. Given a folder, it finds all processed recordings in it and its subfolders, collects metadata and presents it in a table form. It further enables entering of the experimental protocol, and additional data, which are then also searchable. It also provides a link to an automaticaly generated notebook for a brief glimpse into actual results of an experiment. See demo on youtube (https://youtu.be/tj4TjL_PJ1Q).

Installation
============

Requirements
------------
The module is intended to be used within a JupyterHub environment.
It is possible to run in python without a jupyter instance, but features like our apps may not be working in this case.
If you do not want to stress about installing a jupyter environment around it, just use our docker image or use it `bundled in JupyterHub <https://github.com/Hannnsen/Physio_Ca_JupyterHub/>`_.
A java environment is mandatory due to usage of bioformats to convert image data.

- JVM
- Python >=3.8, <3.11
- Modules (and their dependencies) listed in pyproject.toml

PyPI
----
We offer a pre-built version of our software package.
If your system meets the requirements, you can install it easily with a package manager (e.g. pip, poetry, conda) of your choice.
In our hands, poetry was the easiest deployment tool, so we would install it like this into its standard environment:

.. code-block:: sh
    
    poetry install islets


Docker
------
A dockerfile is included in the root of the framework. It contains everything to run python code in the base environment. It can be built with the following command:

.. code-block:: sh

   docker build -t ctn_server .

If you do not want to build it yourself there is a prebuilt version on docker-hub. It can be pulled simly by:

.. code-block:: sh

   docker pull hannsen/cell-tissue-networks_server:latest

As an example, to run the server with custom data and access it in a shell you can use it like this:

.. code-block:: sh
   
   docker run -it -v /path/to/real/data:/data:rw hannsen/cell-tissue-networks_server:latest /bin/bash


.. |Docker Image Version badge| image:: https://img.shields.io/docker/v/hannsen/cell-tissue-networks_server/latest?logo=docker
   :alt: Docker Image Version (tag latest semver)

.. |Docker Image Size badge| image:: https://img.shields.io/docker/image-size/hannsen/cell-tissue-networks_server/latest
   :alt: Docker Image Size (tag)
   
.. |PyPI Version badge| image:: https://img.shields.io/pypi/v/islets
   :alt: PyPI
   
.. |PyPI Python Version badge| image:: https://img.shields.io/pypi/pyversions/islets
   :alt: PyPI - Python Version
