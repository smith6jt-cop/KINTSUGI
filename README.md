# KINTSUGI: Knowledge Integration with New Technologies for Simplified User-Guided Image processing

<p align="center">
  <img src="/docs/CD8_curtain.gif" alt="Description of the gif" style="float: right; margin-left: 20px;">
    
## Multiplex image processing for challenging datasets with a focus on user integration rather than automation.  This pipeline includes 2D/3D GPU/CPU illumination correction, stitching, deconvolution, extended depth of focus, registration, autofluorescence removal, segmentation, clustering, and spatial analysis.
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
    - [5. Pixel Clustering](#5-pixel-clustering)
    - [6. Cell Clustering](#5-cell-clustering)

### Installation Steps

#### 1. Download miniforge 
&emsp;Download and install environment management software.

&emsp;Download miniforge: [https://github.com/conda-forge/miniforge](https://github.com/conda-forge/miniforge.).

&emsp;Follow installation instructions for your OS.

#### 2. Create mamba environment
&emsp;Launch miniforge as administrator (if possible). You will be in the default “base” environment.

&emsp;Change directory to your user folder: 
```
cd C:\Users\[your user name]
cd /home/[your user name]
```
&emsp;For Linux OS, you may need to enter:
```
source "${HOME}/miniforge3/etc/profile.d/mamba.sh"
source "${HOME}/miniforge3/etc/profile.d/conda.sh"
mamba activate
```

&emsp;To download the code and associated files enter: 
```
git clone https://github.com/smith6jt-cop/KINTSUGI.git
```
&emsp;Change directory to enter the folder just downloaded 
```
cd KINTSUGI
```
&emsp;For Windows OS, create the environment by entering:
```
mamba env create -f environment.yml
```
&emsp;For Linux OS, create the environment by entering:
```
mamba env create -f environment_linux.yml
```
&emsp;The downloading and installation of the packages will take several minutes depending on available computing resources and network speed.

&emsp;Activate the environment by entering:
```
mamba activate KINTSUGI
```
&emsp;It is recommended to use VS Code to run the notebooks. Download and install VS Code [https://code.visualstudio.com/](https://code.visualstudio.com/).


#### 3. Download files
&emsp;Download zip files and extract them to KINTSUGI folder. 

&emsp;&emsp;Java - Information at: [https://www.oracle.com/java/technologies/downloads/#java21](https://www.oracle.com/java/technologies/downloads/#java21). 
  
&emsp;&emsp;Download links:  

&emsp;&emsp;&emsp;[https://download.oracle.com/java/21/latest/jdk-21_windows-x64_bin.zip (sha256)](https://download.oracle.com/java/21/latest/jdk-21_windows-x64_bin.zip)   
&emsp;&emsp;&emsp;[https://download.oracle.com/java/21/latest/jdk-21_macos-aarch64_bin.tar.gz (sha256)](https://download.oracle.com/java/21/latest/jdk-21_linux-x64_bin.tar.gz)   
&emsp;&emsp;&emsp;[https://download.oracle.com/java/21/latest/jdk-21_linux-aarch64_bin.tar.gz (sha256) ](https://download.oracle.com/java/21/latest/jdk-21_linux-aarch64_bin.tar.gz)  

&emsp;&emsp;Maven - Information at: [https://maven.apache.org/download.cgi.](https://maven.apache.org/download.cgi) 

&emsp;&emsp;Download links:   

&emsp;&emsp;&emsp;[apache-maven-3.9.9-bin.zip ](https://dlcdn.apache.org/maven/maven-3/3.9.9/binaries/apache-maven-3.9.9-bin.zip)  
&emsp;&emsp;&emsp;[apache-maven-3.9.9-bin.tar.gz ](https://dlcdn.apache.org/maven/maven-3/3.9.9/binaries/apache-maven-3.9.9-bin.tar.gz)  

&emsp;&emsp;PyVips (LibVips) (for VALIS Registration only).   

&emsp;&emsp;&emsp;Windows download link: [vips-dev-w64-all-8.16.0.zip ](https://github.com/libvips/build-win64-mxe/releases/download/v8.16.0/vips-dev-w64-all-8.16.0.zip)  
&emsp;&emsp;&emsp;Additional install instructions for Linux: [https://github.com/libvips/pyvips](https://github.com/libvips/pyvips)    

&emsp;&emsp;FIJI/ImageJ: [https://imagej.net/software/fiji/downloads](https://imagej.net/software/fiji/downloads)

&emsp;&emsp;&emsp;Install the "Fiji.app" folder to your user folder.  
&emsp;&emsp;&emsp;Follow clij2 installation: [https://clij.github.io/clij2-docs/installationInFiji](https://clij.github.io/clij2-docs/installationInFiji)  


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

[5. Pixel Clustering](notebooks/5_Cluster_Pixels.ipynb)
  For self-organizing map application to pixels.

[6. Cell Clustering](notebooks/6_Cluster_Cells.ipynbb)
  For self-organizing map application to pixel clusters and segmetation features.

<div>

## Acknowledgements

Shoulders of giants we stand on:

&emsp;For multiplex histology/imaging: Nolan lab  
&emsp;For illumination correction: Peng lab  
&emsp;For stitching: MIST, Fukai's m2stitch  
&emsp;For deconvolution: Becker's LsDeconv  
&emsp;For EDoF: Forster et al., Clij2  
&emsp;For registration: Gatenbee's VALIS  
&emsp;For segmentation: VanValen lab  
&emsp;For clustering: Angelo lab  
&emsp;For general processing: ImageJ/FIJI, pyImageJ, Haase's stackview  
&emsp;Python, Jupyter, conda/mamba, Java, Maven, Linux, Windows, UF, NIH, HubMap, the power grid, the internet, the earth, gravity, the sun, oxygen, caffeine, love, neurons, neurotransmitters, water, and the unkown.
