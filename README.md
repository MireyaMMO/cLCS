# cLCS
## Description
Climatological LCS based on the code developed by Rodrigo Duran for [Matlab](https://bitbucket.org/rodu/clcss/src/master/)

This repository contains the scripts needed to deploy a series of particle releases using OpenDrift and calculate the  associated climatological Lagrangian Coherent Structures.


## Requirements
To use this package OpenDrift is necessary. To facilitate and make the installation faster please use [mamba](https://github.com/conda-forge/miniforge). 
As shown in the link to install on Unix-like platforms (Mac OS & Linux) use:
```
$ curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh

$ bash Miniforge3-$(uname)-$(uname -m).sh
```

Once mamba is installed we can proceed.

This research was part of the [Moana Project](www.moanaproject.org), therefore the reader used for OpenDrift is especific to a reader developed by [Simon Weppe](https://github.com/simonweppe). However, unless working with the same dataset or a data related to the Moana project I invite you to download the repository directly from OpenDrift

OpenDrift original repo
```
$ git clone https://github.com/OpenDrift/opendrift.git
```

Branch developed by Simon Weppe
```
$ git clone https://github.com/simonweppe/opendrift.git
```

```
$ cd ../
$ git clone https://github.com/MireyaMMO/cLCS.git 
$ cd cLCS
$ mamba env create --name cLCS --file=environment.yml
$ pip install --no-deps -e .
$ cd ../opendrift
$ pip install --no-deps -e .

```

## Contents
- cLCS
  - make_cLCS.py 
  - mean_C.py: meanC computes monthly averages for each month
  - plotting.py

- Examples
  - jupyter_notebooks with examples 

## References:
[Duran, R., Beron-Vera, F.J. & Olascoaga, M.J. Extracting quasi-steady Lagrangian transport patterns from the ocean circulation: An application to the Gulf of Mexico. Sci Rep 8, 5218 (2018). https://doi.org/10.1038/s41598-018-23121-y](https://www.nature.com/articles/s41598-018-23121-y)

[Duran, R., F. J. Beron-Vera and M. J. Olascoaga (2019). Climatologial Lagrangian Coherent Structures code. DOI: 10.18141/1558781](https://bitbucket.org/rodu/clcss/src/master/)

[Dagestad, K.-F., Rohrs, J., Breivik, O., and Ådlandsvik, B. (2018). OpenDrift v1.0: a generic framework for trajectory modelling. Geoscientific Model Development, 11 (4), 1405–1420.](https://github.com/OpenDrift/opendrift)


Montano Orozco, M. M. (2023). Hydrodynamics and coastal dispersion from the Bay of Plenty, Aotearoa, New Zealand: A 25-year numerical modelling perspective (Thesis, Doctor of Philosophy). University of Otago. Retrieved from http://hdl.handle.net/10523/16243. 

[For a better formatted version of the thesis (working embedded links) please click here](https://drive.google.com/file/d/1WMgq2lu7K5MjGTy6O5YpoKDkONDclkHo/view?usp=sharing)