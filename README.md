# cLCS
## Description
Climatological LCS based on the code developed by Rodrigo Duran for Matlab

This repository contains the scripts needed to deploy a series of particle releases and calculate the climatological Lagrangian Coherent Structures.

## Requirements
Please first install OpenDrift following the instructions
```
$ git clone https://github.com/OpenDrift/opendrift.git
$ cd opendrift/
$ conda env create --name OpenDrift --file=environment.yml
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

[OpenDrift](https://github.com/OpenDrift/opendrift)