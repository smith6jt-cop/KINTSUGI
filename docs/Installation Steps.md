# Installation Steps

## 1. Download miniforge 
&emsp;Download and install environment management software.

&emsp;Download miniforge: [https://github.com/conda-forge/miniforge](https://github.com/conda-forge/miniforge.).

&emsp;Follow installation instructions for your OS.

## 2. Create mamba environment
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


## 3. Download files
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


## 4. Copy/move raw image data  
&emsp;Create a folder in the KINTSUGI folder called “data”.  

&emsp;If downloading test data use this link: [https://uflorida-my.sharepoint.com/:f:/g/personal/smith6jt_ufl_edu1/Er5ui-wFA6BNnmgj9N1hPAsBYQaiKfSQa2do_lUMhQdaGg?e=5Uny95](https://uflorida-my.sharepoint.com/:f:/g/personal/smith6jt_ufl_edu1/Er5ui-wFA6BNnmgj9N1hPAsB_Z8EwL7jkfekJwrWEfVRbw?e=oxaxMH)  

&emsp;Move all image data to [your user folder]\KINTSUGI\data.  


<div>