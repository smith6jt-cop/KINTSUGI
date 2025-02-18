# KINTSUGI: Knowledge Integration with New Technologies for Simplified User-Guided Image processing

<p align="center">
  <img src="/docs/CD31_curtain.gif" alt="CD31 Autofluorescence Removal" style="float: right; margin-left: 20px;">
    
## Multiplex image processing for challenging datasets with a focus on user integration rather than automation.  This pipeline includes 2D/3D GPU/CPU illumination correction, stitching, deconvolution, extended depth of focus, registration, autofluorescence removal, segmentation, clustering, and spatial analysis.

To view notebooks in browser see: [https://smith6jt-cop.github.io/KINTSUGI-docs/Intro.html](https://smith6jt-cop.github.io/KINTSUGI-docs/Intro.html)

</p>

<div>
  
## Table of Contents

  - [Installation Steps](#installation-steps)
    - [1. Download miniforge](#1-download-miniforge)
    - [2. Create mamba environment](#2-create-mamba-environment)
    - [3. Download files](#3-download-files)
    - [4. Copy/move raw image data](#4-copy/move-raw-image-data)
  - [Notebooks](#notebooks)
    - [1. Parameter tuning/testing](#1-parameter-tuning/testing)
    - [2. Batch processing](#2-batch-processing)
    - [3. Signal Isolation](#3-signal-isolation)
    - [4. Segmentation](#4-segmentation)

### Installation Steps

#### 1. Download miniforge 
&emsp;Download and install conda-based environment management software.  If you already have, skip to step 2.

&emsp;Download miniforge: [https://github.com/conda-forge/miniforge](https://github.com/conda-forge/miniforge.).

&emsp;Follow installation instructions for your OS.

#### 2. Create mamba environment
&emsp;Launch miniforge as administrator (if possible). You will be in the default “base” environment.

&emsp;Change directory to your user folder: 
```
cd C:\Users\[your user name]
```
&emsp;To download the code and associated files enter: 
```
git clone https://github.com/smith6jt-cop/KINTSUGI.git
```
&emsp;Change directory to enter the folder just downloaded 
```
cd KINTSUGI
```
&emsp;Create the environment by entering:
```
mamba env create -f environment.yml
```
```
&emsp;Activate the environment by entering:
```
mamba activate KINTSUGI
```
&emsp;It is recommended to use VS Code to run the notebooks. Download and install VS Code [https://code.visualstudio.com/](https://code.visualstudio.com/).


#### 3. Download dependency files
&emsp;Download zip files and extract them to KINTSUGI folder. 

#### 4. Copy/move raw image data  
&emsp;Create a folder in the KINTSUGI folder called “data”.  

&emsp;If downloading test data use this link: [https://uflorida-my.sharepoint.com/:f:/g/personal/smith6jt_ufl_edu1/Er5ui-wFA6BNnmgj9N1hPAsBYQaiKfSQa2do_lUMhQdaGg?e=5Uny95](https://uflorida-my.sharepoint.com/:f:/g/personal/smith6jt_ufl_edu1/Er5ui-wFA6BNnmgj9N1hPAsB_Z8EwL7jkfekJwrWEfVRbw?e=oxaxMH)  

&emsp;Move all image data to [your user folder]\KINTSUGI\data.  


<div>


## Notebooks
[1. Parameter tuning/testing](notebooks/1_Single_Channel_Eval.ipynb)
  For testing illumination correction, stitching, deconvolution, and EDoF.

[2. Batch processing](notebooks/2_Cycle_Processing.ipynb)
  For batch processing illumination correction, stitching, deconvolution, EDoF, and registration.

[3. Signal Isolation](notebooks/3_Signal_Isolation.ipynb)
  For autofluorescence subtraction, filtering, and final processing to isolate signal.

[4. Segmentation](notebooks/4_Segmentation.ipynb)
  For Mesmer segmentation and feature extraction.

<div>

