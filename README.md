# KINTSUGI: Knowledge Integration with New Technologies for Simplified User-Guided Image processing

<p align="center">
  <img src="/docs/CD3e.gif" alt="CD31 Autofluorescence Removal" style="float: right; margin-left: 20px;">
    
## Multiplex image processing for challenging datasets with a focus on user integration rather than automation.  This pipeline includes 2D/3D GPU/CPU illumination correction, stitching, deconvolution, extended depth of focus, registration, autofluorescence removal, segmentation, clustering, and spatial analysis.

Citation Information:

Smith, J. A. et al. Protocol for processing and analyzing multiplexed images improves lymphatic cell identification and spatial architecture in human tissue. STAR Protocols 6, 103976 (2025).

</p>

[![DOI](https://zenodo.org/badge/794118146.svg)](https://doi.org/10.5281/zenodo.14984518)

<div>
  
## Table of Contents

  - [Installation Steps](#installation-steps)
    - [1. Download miniconda](#1-download-miniconda)
    - [2. Create environment](#2-create-environment)
    - [3. Download dependency files](#3-download-dependency-files)
    - [4. Copy or move raw image data](#4-copy-or-move-raw-image-data)
    - [5. Setup VS Code](#5-setup-vs-code)
  - [Notebooks](#notebooks)
    - [1. Parameter tuning and testing](#1-parameter-tuning-and-testing)
    - [2. Batch processing](#2-batch-processing)
    - [3. Signal Isolation](#3-signal-isolation)
    - [4. Segmentation](#4-segmentation)

### Installation Steps

#### 1. Download miniconda 
*  Download and install conda-based environment management software.  If you already have, skip to step 2.

*  Download miniconda: [https://www.anaconda.com/download/success#miniconda](https://www.anaconda.com/download/success#miniconda).

*  Follow installation instructions.

#### 2. Create environment
*  Open a miniconda "Anaconda Prompt” terminal.   You will be in the default “base” environment. Update the base environment and to set the default dependency solver to libmamba by copying and entering the following: 
```
conda update -n base conda 
conda install -n base conda-libmamba-solver 
conda config --set solver libmamba
```
*  Change directory to your user folder: 
```
cd C:\Users\[your user name]
```
*  Install “git” to the base environment: 
```
conda install git -y 
```
*  To download the code and associated files enter: 
```
git clone https://github.com/smith6jt-cop/KINTSUGI.git
```
*  Change directory to enter the folder just downloaded 
```
cd KINTSUGI
```
*  Create the environment by entering:
```
conda env create -f env.yml
```
*  Close the Anaconda PowerShell Prompt 

#### 3. Download dependency files

*  There are necessary files that are too large to host on GitHub. To download the extra dependency zip file, use this link: [https://zenodo.org/records/14969214](https://zenodo.org/records/14969214?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6ImVkNGE2ZWQxLWNjNGItNDgzNi04MTBmLWI0ZmU2MWM1N2Q1MSIsImRhdGEiOnt9LCJyYW5kb20iOiJjYmUwY2U2YTg1ODc1YmQ4MGE2NDk4NjI4ZDQ1ZTcwYSJ9.uMgV9EeOeGi6MsaPqrVgasA1oNoDI7SWFtKU6OK_RZ_BdOpHsq-4XHB-jfQA9aV3tVWbp23cD2XqL4B3VSEiEw) and extract contents.  Extract each zipped file to the KINTSUGI folder.  Alternatively, you may download, install, and configure maven3.9.9, java-jdk21, MatlabRuntime 2024a, PyVips-dev-8.16, and FIJI with the Clij2 plugin.

*  Download zip files and extract them to KINTSUGI folder. 

#### 4. Copy or move raw image data  
*  Create a folder in the KINTSUGI folder called “data”.  

*  To download test data, use the KINTSUGI Zenodo community: [https://zenodo.org/communities/kintsugi](https://zenodo.org/communities/kintsugi/records?q=&l=list&p=2&s=10&sort=newest)

*  Results of processing the test dataset can be found at: [https://app.globus.org]( https://app.globus.org/file-manager?origin_id=10f408d9-f5ee-11ef-bf21-0affeb6b961d&origin_path=%2F)

*  Move all image data to [your user folder]\KINTSUGI\data.  

#### 5. Setup VS Code
*  It is recommended to use VS Code to run the notebooks. Download and install VS Code [https://code.visualstudio.com/](https://code.visualstudio.com/).

*  Create a GitHub account if you do not have one.  This is necessary for using the AI assistant CoPilot in VSCode. 

*  Launch VS Code and sign in to GitHub. 

*  Close VS Code.  

*  Launch Anaconda PowerShell Prompt 

*  Change the directory: 
```
cd C:\Users\[username]\KINTSUGI
```
*  Activate the environment by entering: 
```
conda activate KINTSUGI
```
*  Launch VS Code by entering “code .” (the word ‘code’ followed by a space and a period) from your conda terminal in the activated KINTSUGI environment. 

*  Launching VS Code from the activated environment in the Anaconda Prompt terminal is the only way KINTSUGI should be initiated.  This will ensure all package functions are available.

*  If prompted to trust the authors of the files in the folder, check the box and click “Yes, I trust the authors.” 

*  If prompted to Download or configure “Git” click “Do Not Show Again”. 

*  If prompted to “install the recommended extensions ...” click “Install”. 

<div>

### Notebooks
#### 1. Parameter tuning and testing
*  For testing illumination correction, stitching, deconvolution, and EDoF. [notebooks/1_Single_Channel_Eval.ipynb](notebooks/1_Single_Channel_Eval.ipynb)

#### 2. Batch processing
*  For batch processing illumination correction, stitching, deconvolution, EDoF, and registration. [notebooks/2_Cycle_Processing.ipynb](notebooks/2_Cycle_Processing.ipynb)

#### 3. Signal Isolation
*  For autofluorescence subtraction, filtering, and final processing to isolate signal.[notebooks/3_Signal_Isolation.ipynb](notebooks/3_Signal_Isolation.ipynb)

#### 4. Segmentation
*  For Instanseg segmentation, feature extraction, and spatial analysis.[notebooks/4_Segmentation.ipynb](notebooks/4_Segmentation.ipynb)

<div>

