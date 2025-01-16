# KINTSUGI: Knowledge Integration with New Technologies for Simplified User-Guided Image processing

Multiplex image processing for challenging datasets with a focus on user integration rather than automation.

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
Download and install environment management software.

Download miniforge: https://github.com/conda-forge/miniforge.

Follow installation instructions for your OS.

#### 2. Create mamba environment
Launch miniforge as administrator (if possible). 

You will be in the default “base” environment.

If Windows OS, change directory to your user folder by entering without the quotes: 
```sh
cd C:\Users\[your user name]\
```
For Linux OS, open a terminal and enter:
```
source "${HOME}/miniforge3/etc/profile.d/mamba.sh"
source "${HOME}/miniforge3/etc/profile.d/conda.sh"
mamba activate
```

To download the code and associated files enter: 
```
git clone https://github.com/smith6jt-cop/KINTSUGI.git
```
Change directory to enter the folder just downloaded 
```
cd KINTSUGI
```
For Windows OS, create the environment by entering:
```
mamba env create -f environment.yml
```
For Linux OS, create the environment by entering:
```
mamba env create -f environment_linux.yml
```
The downloading and installation of the packages will take several minutes depending on available computing resources and network speed.

Activate the environment by entering:
```
mamba activate KINTSUGI
```
It is recommended to use VS Code to run the notebooks. Download and install VS Code [https://code.visualstudio.com/](https://code.visualstudio.com/).

#### 3. Download files
Download zip files and extract them to KINTSUGI folder 

  Java - Information at: [https://www.oracle.com/java/technologies/downloads/#java21](https://www.oracle.com/java/technologies/downloads/#java21). Download links:  

  [https://download.oracle.com/java/21/latest/jdk-21_windows-x64_bin.zip (sha256)](https://download.oracle.com/java/21/latest/jdk-21_windows-x64_bin.zip) 

  [https://download.oracle.com/java/21/latest/jdk-21_macos-aarch64_bin.tar.gz (sha256)](https://download.oracle.com/java/21/latest/jdk-21_linux-x64_bin.tar.gz) 

  [https://download.oracle.com/java/21/latest/jdk-21_linux-aarch64_bin.tar.gz (sha256) ](https://download.oracle.com/java/21/latest/jdk-21_linux-aarch64_bin.tar.gz)

Maven - Information at: [https://maven.apache.org/download.cgi.](https://maven.apache.org/download.cgi) Download links:   

  [apache-maven-3.9.9-bin.zip ](https://dlcdn.apache.org/maven/maven-3/3.9.9/binaries/apache-maven-3.9.9-bin.zip)

  [apache-maven-3.9.9-bin.tar.gz ](https://dlcdn.apache.org/maven/maven-3/3.9.9/binaries/apache-maven-3.9.9-bin.tar.gz)

PyVips (LibVips) (for VALIS Registration only).   

  Windows download link: [vips-dev-w64-all-8.16.0.zip ](https://github.com/libvips/build-win64-mxe/releases/download/v8.16.0/vips-dev-w64-all-8.16.0.zip)

  Additional install instructions for Linux: [https://github.com/libvips/pyvips](https://github.com/libvips/pyvips)  

FIJI/ImageJ: [https://imagej.net/software/fiji/downloads](https://imagej.net/software/fiji/downloads)

  Install the "Fiji.app" folder to your user folder.

  Follow clij2 installation: [https://clij.github.io/clij2-docs/installationInFiji](https://clij.github.io/clij2-docs/installationInFiji)

#### 4. Copy/move raw image data
Create a folder in the KINTSUGI folder called “data”.

If downloading test data use this link: [https://uflorida-my.sharepoint.com/:f:/g/personal/smith6jt_ufl_edu1/Er5ui-wFA6BNnmgj9N1hPAsBYQaiKfSQa2do_lUMhQdaGg?e=5Uny95](https://uflorida-my.sharepoint.com/:f:/g/personal/smith6jt_ufl_edu1/Er5ui-wFA6BNnmgj9N1hPAsB_Z8EwL7jkfekJwrWEfVRbw?e=oxaxMH)

Move all image data to [your user folder]\KINTSUGI\data.


<div>


### Notebooks
[1. Parameter tuning/testing](notebooks/1_Single_Channel_Eval.ipynb) 

[2. Batch processing](notebooks/2_Cycle_Processing.ipynb) 

[3. Signal Isolation](notebooks/3_Signal_Isolation.ipynb)

[4. Segmentation](notebooks/4_Segmentation.ipynb)

[5. Pixel Clustering](notebooks/5_Cluster_Pixels.ipynb)

[6. Cell Clustering](notebooks/6_Cluster_Cells.ipynbb)
