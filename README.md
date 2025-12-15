# Introduction

_dumux-rosi_ provides Python interfaces for solving the Richards equation as well as general advection–diffusion–reaction systems in plant–soil interaction modelling. It leverages the nonlinear finite-volume solver of [DuMu<sup>x</sup>](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux) and is designed to integrate with [CPlantBox](https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox).

# Installation

## Using a Pyton script 

The installation scipt is located within the [CPlantBox](https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox) repository.
This script will install the full setup:  CPlantBox,  _dumux-rosi_, and DuMu<sup>x</sup>]. 
Just download and run the Python file "installDumuxRosi_Ubuntu.py" (which is based on the DuMu<sup>x</sup> installation file).
```bash
sudo apt-get update
sudo apt-get upgrade
[ ! -d 'cpbenv' ] && python3 -m venv cpbenv &&  source cpbenv/bin/activate ||  source cpbenv/bin/activate
wget https://raw.githubusercontent.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox/master/installDumuxRosi_Ubuntu.py
python3 installDumuxRosi_Ubuntu.py
```
Finally, run the installation script from the dumux-rosi directory:
```bash
cd dumux-rosi
./install_modules.sh
```
This will install dumux-rosi into your Python site-packages, making it available for import.

## by hand

* Use the DuMu<sup>x</sup> [Installation Notes](https://dumux.org/docs/doxygen/master/installation.html) to set up DuMu<sup>x</sup>
* Install [CPlantBox](https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox) into the DuMu<sup>x</sup>  base folder (where dune-common/, dumux/, etc. are located) according to the CPlantBox installation instruction (README.md)
* Finally, clone the _dumux-rosi_ repository into the DuMu<sup>x</sup> common base folder and use DUNE build system to build the repository, and run install_modules.sh in the dumux-rosi folder to make it ready for use in Python.



# Further 

Please read the manual (Manual.pdf on this repository) 
- for installation guidelines
- getting started with examples
- folder structure is explained in the dumux-rosi wiki https://github.com/Plant-Root-Soil-Interactions-Modelling/dumux-rosi/wiki/Folder-structure

<img src="Logo_long_white.png" alt="drawing" width="200"/>
