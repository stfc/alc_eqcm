## About the code
**ALC_EQCM** is a software for post-processing Electrochemical Quartz Crystal Microbalance (EQCM) data. The implemented functionalities allow quantitative characterization of electrochemical processes and interpretation of stoichiometric changes. This information is used to automatically generate atomistic models compatible with EQCM data, together with appropriate settings for atomistic level simulations of the derived models. These software capabilities constitute an alternative tool that aims to bridge the gap between experimental and computational research of complex electrochemical reactions.

The development of this code started in April 2020 at the [Ada Lovelace Centre](https://adalovelacecentre.ac.uk/) (ALC) of the [Science and Technology Facilities Council](https://stfc.ukri.org/) (STFC). **ALC_EQCM** is a serial code written in modern Fortran according to the 2008 standards. Its structure for development and maintenance follows the Continuous Integration (CI) practice and it is integrated within the GitLab DevOps of the STFC.

Together with the code, we also provide a manual intended to guide the user with the multiple functionalities and help with the initial challenges for compilation and execution. We remind the user that the manual is a living document, which is subject to continuous modifications and updates. For this reason, we strongly advise the user to refer to the most updated version of the manual.

## License
ALC_EQCM redistribution is under the BSD 3-Clause License. Please refer to file [LICENSE](./LICENSE) in the root directory for further details.

## Structure of files and folders
ALC_EQCM contains the following set of files and folders (in italic-bold):

* [***CI-tests***](./CI-tests): contains all the tests in .tar format for testing purposes. There is also a file called *README.txt* with a brief description for each test.
* [***cmake***](./cmake): contains the specification for the compilation flags depending on the Fortran compiler, including options for debugging.
* [cmake_building.md](./cmake_building.md): details the steps to build, compile and run tests using the CMake platform.
* [CMakeList.txt](./CMakeList.txt): sets the framework for code building and testing with CMake.
* [LICENSE](./LICENSE): conditions for ALC_EQCM license.
* [***manual***](./manual): folder with the user's manual.
* README.md: this file .
* [***scripts***](./scripts): contains scripts for data processing and transformation of atomistic structures with the *.cif* format
* [***source***](./source): contains the source code. Files have the *.F90* extension
* [***tools***](./tools): includes all shell files for building, compiling and testing the code automatically.
* [***tutorial***](./tutorials): contains the input files for the tutorial examples, which are explained in the section 6 of the manual

## Contributors
 * Ivan Scivetti
 * Gilberto Teobaldi 

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

### Depedencies
The user must have access to the following software:  

* GNU-Fortran (5.4.0) or Intel-Fortran (16.0.1)
* Cmake (3.1)  
* Make (3.82)  
* git (2.7.4)  
The following two softwares are only needed when working with atomic structures in *.cif* format:
* Python (2.7.12)
* Atomistic Simulation Environment-ASE (2.1)

Information in parenthesis indicates the minimum version tested during the development of the code. The specification for the minimum versions is not fully rigorous but indicative, as there could be combinations of other minimum versions that still work. Our tests indicate that versions of Intel compiler older than 16.0.1 exhibit problems and should be avoided.

### Building and testing the code with CMake
Details can be found in file [cmake_building.md](./cmake_building.md)

