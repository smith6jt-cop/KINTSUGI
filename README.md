# KINTSUGI: Knowledge Integration with New Technologies for Simplified User-Guided Image processing

<p align="center">
  <img src="/docs/CD31_curtain.gif" alt="CD31 Autofluorescence Removal" style="float: right; margin-left: 20px;">
    
## Multiplex image processing for challenging datasets with a focus on user integration rather than automation.  This pipeline includes 2D/3D GPU/CPU illumination correction, stitching, deconvolution, extended depth of focus, registration, autofluorescence removal, segmentation, clustering, and spatial analysis.


</p>

<div>
  
## Table of Contents

  - [Installation Steps](#installation-steps)
    - [1. Download miniconda](#1-download-miniconda)
    - [2. Create environment](#2-create-environment)
    - [3. Download dependency files](#3-download-dependency-files)
    - [4. Copy/move raw image data](#4-copy/move-raw-image-data)
    - [5. Setup VS Code](#5-setup-vs-code)
  - [Notebooks](#notebooks)
    - [1. Parameter tuning/testing](#1-parameter-tuning/testing)
    - [2. Batch processing](#2-batch-processing)
    - [3. Signal Isolation](#3-signal-isolation)
    - [4. Segmentation](#4-segmentation)

### Installation Steps

#### 1. Download miniconda 
&emsp;Download and install conda-based environment management software.  If you already have, skip to step 2.

&emsp;Download miniconda: [https://www.anaconda.com/download/success#miniconda](https://www.anaconda.com/download/success#miniconda).

&emsp;Follow installation instructions.

#### 2. Create environment
&emsp;Open a miniconda "Anaconda Prompt” terminal.   You will be in the default “base” environment. Update the base environment and to set the default dependency solver to libmamba by copying and entering the following: 
```
conda update -n base conda 
conda install -n base conda-libmamba-solver 
conda config --set solver libmamba
```
&emsp;Change directory to your user folder: 
```
cd C:\Users\[your user name]
```
&emsp;Install “git” to the base environment: 
```
conda install git -y 
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
conda env create -f env.yml
```

&emsp;Close the Anaconda PowerShell Prompt 

#### 3. Download dependency files

&emsp;There are necessary files that are too large to host on GitHub. To download the extra dependency zip file, use this link: [https://drive.google.com/file/d/1_CoQ2o4iqc1HT9AlrlmrawKEt5A9dWJx/view?usp=sharing](https://drive.google.com/file/d/1_CoQ2o4iqc1HT9AlrlmrawKEt5A9dWJx/view?usp=sharing) and extract contents.  Extract each zipped file to the KINTSUGI folder.  Alternatively, you may download, install, and configure maven3.9.9, java-jdk21, MatlabRuntime 2024a, PyVips-dev-8.16, and FIJI with the Clij2 plugin.

&emsp;Download zip files and extract them to KINTSUGI folder. 

#### 4. Copy/move raw image data  
&emsp;Create a folder in the KINTSUGI folder called “data”.  

&emsp;To download test data, use Globus with your institution or Google/ORCID/Github account: [https://app.globus.org/file-manager?origin_id=dce1f3d9-f067-11ef-a905-0207be7ee3a1&origin_path=%2F](https://app.globus.org/file-manager?origin_id=dce1f3d9-f067-11ef-a905-0207be7ee3a1&origin_path=%2F) Alternatively download from the KINTSUGI Zenodo community: [https://zenodo.org/communities/kintsugi/records?q=&l=list&p=2&s=10&sort=newest](https://zenodo.org/communities/kintsugi/records?q=&l=list&p=2&s=10&sort=newest)

Results of processing the test dataset can be found at: [ https://app.globus.org/file-manager?origin_id=10f408d9-f5ee-11ef-bf21-0affeb6b961d&origin_path=%2F]( https://app.globus.org/file-manager?origin_id=10f408d9-f5ee-11ef-bf21-0affeb6b961d&origin_path=%2F)

&emsp;Move all image data to [your user folder]\KINTSUGI\data.  

#### 5. Setup VS Code
&emsp;It is recommended to use VS Code to run the notebooks. Download and install VS Code [https://code.visualstudio.com/](https://code.visualstudio.com/).

&emsp;Create a GitHub account if you do not have one.  This is necessary for using the AI assistant CoPilot in VSCode. 

&emsp;Launch VS Code and sign in to GitHub. 

&emsp;Close VS Code.  

&emsp;Launch Anaconda PowerShell Prompt 

&emsp;Change the directory: 
```
cd C:\Users\[username]\ KINTSUGI
```
&emsp;Activate the environment by entering: 
```
conda activate KINTSUGI
```
&emsp;Launch VS Code by entering “code .” (the word ‘code’ followed by a space and a period) from your conda terminal in the activated KINTSUGI environment. 

&emsp;Launching VS Code from the activated environment in the Anaconda Prompt terminal is the only way KINTSUGI should be initiated.  This will ensure all package functions are available.

&emsp;If prompted to trust the authors of the files in the folder, check the box and click “Yes, I trust the authors.” 

&emsp;If prompted to Download or configure “Git” click “Do Not Show Again”. 

&emsp;If prompted to “install the recommended extensions ...” click “Install”. 
<div>


## Notebooks
[1. Parameter tuning/testing](notebooks/1_Single_Channel_Eval.ipynb)
  For testing illumination correction, stitching, deconvolution, and EDoF.

[2. Batch processing](notebooks/2_Cycle_Processing.ipynb)
  For batch processing illumination correction, stitching, deconvolution, EDoF, and registration.

[3. Signal Isolation](notebooks/3_Signal_Isolation.ipynb)
  For autofluorescence subtraction, filtering, and final processing to isolate signal.

[4. Segmentation](notebooks/4_Segmentation.ipynb)
  For Instanseg segmentation, feature extraction, and spatial analysis.

<div>

