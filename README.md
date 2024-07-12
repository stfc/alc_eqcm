## Project status
Funding support for full-time development is not longer available. Code development now relies on the generosity of those individuals who want to contribute. If you are interested, please contact Ivan Scivetti (ivan.scivetti@stfc.ac.uk) or Gilberto Teobaldi (gilberto.teobaldi@stfc.ac.uk).  
Alternatively, you can contact the Ada Lovelace Centre (ALC) via the [website](https://adalovelacecentre.ac.uk/contact-us/) or via email (alc@stfc.ac.uk).  

## Description
**ALC_EQCM** is a software for post-processing Electrochemical Quartz Crystal Microbalance (EQCM) data. The implemented functionalities allow quantitative characterization of electrochemical processes and interpretation of stoichiometric changes. This information is used to automatically generate atomistic models compatible with EQCM data, together with appropriate settings for atomistic level simulations of the derived models. These software capabilities constitute an alternative tool that aims to bridge the gap between experimental and computational research of complex electrochemical reactions.

The development of this code started in April 2020 at the [ALC](https://adalovelacecentre.ac.uk/) of the [Science and Technology Facilities Council](https://stfc.ukri.org/) (STFC). **ALC_EQCM** is a serial code written in modern Fortran according to the 2008 standards. Its structure for development and maintenance follows the Continuous Integration (CI) practice and it is integrated within the GitLab DevOps of the STFC.

Together with the code, we also provide a manual intended to guide the user with the multiple functionalities and help with the initial challenges for compilation and execution. We remind the user that the manual is a living document, which is subject to continuous modifications and updates. For this reason, we strongly advise the user to refer to the most updated version of the manual.

## Disclaimer
Ada Lovelace Centre does not fully guarantee the code is free of errors and assumes no legal responsibility for any incorrect outcome or loss of data.

## Citing ALC_EQCM
Please cite the following work in publications making use of ALC_EQCM:

1) Quantitative Resolution of Complex Stoichiometric Changes During Electrochemical Cycling by Density Functional Theory Assisted, Electrochemical Quartz Crystal Microbalance. T-H. Wu; I. Scivetti; J-C. Chen; J-A. Wang; G. Teobaldi, C-C Hu; L.J. Hardwick. ACS Appl. Energy Mater. 3, 4, 3347–3357 (2020), https://doi.org/10.1021/acsaem.9b02386

2) ALC_EQCM: Automated stoichiometric resolution in electrochemistry through Density Functional Theory aided, Electrochemical Quartz Crystal Microbalance. I. Scivetti and G. Teobaldi. Computational Materials Science 218, 111968, (2023), https://doi.org/10.1016/j.commatsci.2022.111968

Both references are provided in bibtex format within the [***biblio-references***](./biblio-references) folder in the root directory.

## Structure of files and folders
ALC_EQCM contains the following set of files and folders (in italic-bold):

* [***CI-tests***](./CI-tests): contains all the tests in .tar format for testing purposes. There is also a file called *README.txt* with a brief description for each test.
* [***biblio-references***](./biblio-references): includes the bibliographic references of the ALC_EQCM code in BibTex format. 
* [***cmake***](./cmake): contains the specification for the compilation flags depending on the Fortran compiler, including options for debugging.
* [***manual***](./manual): folder with the user's manual.
* [***scripts***](./scripts): contains scripts for data processing and transformation of atomistic structures with the *.cif* format
* [***source***](./source): contains the source code. Files have the *.F90* extension
* [***tools***](./tools): includes all shell files for building, compiling and testing the code automatically.
* [***tutorial***](./tutorials): contains the input files for the tutorial examples, which are explained in the section 6 of the manual
* [.gitignore](./.gitignore): instructs Git which file to ignore for development and integration.
* [CMakeList.txt](./CMakeList.txt): sets the framework for code building and testing with CMake.
* Jenkinsfile: file with specifications to build and run the testing infrastructure.
* LICENSE: specification of the BSD 3-Clause License under which ALC_EQCM is registered.
* README.md: this file.
* [cmake_building.md](./cmake_building.md): details the steps to build, compile and run tests using the CMake platform.

## Contributors
 * Ivan Scivetti (original author)
 * Gilberto Teobaldi (project management and scientific support)

## Getting started  
### Depedencies
The user must have access to the following software:  

* GNU-Fortran (5.4.0) or Intel-Fortran (16.0.1)
* Cmake (3.1)  
* Make (3.82)  
* git (2.7.4)  
The following two softwares are only needed when working with atomic structures in *.cif* format:
* Python (3.8.10)
* Atomistic Simulation Environment-ASE (ase-3.23.0b1)

Information in parenthesis indicates the minimum version tested during the development of the code. The specification for the minimum versions is not fully rigorous but indicative, as there could be combinations of other minimum versions that still work. Our tests indicate that versions of Intel compiler older than 16.0.1 exhibit problems and should be avoided.

## Getting started
The user with account *"username"* can clone the code locally (in machine *"wherever"*) by executing the following command with the SSH protocol
```sh
username@wherever:/home/username/codes$ git clone git@github.com:stfc/alc_eqcm.git
```
Instead, if the user wants to use the HTTPS protocol it must execute
```sh
username@wherever:/home/username/codes$ git clone https://github.com/stfc/alc_eqcm.git
```
Both ways generate the ***alc_eqcm*** folder as the root directory. Alternatively, the code can be downloaded from any of the available assets.


### Building and testing the code with CMake
Details can be found in file [cmake_building.md](./cmake_building.md)
