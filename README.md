# Abaqus Non-linear Unit Cell and Equivalent Single Layer Scripts 

This repository contains the Python and Fortran scripts, as well as Abaqus input files used in [1]. The code is newer version of what was developed in my master's thesis [2] and should be regarded as the superior version. 

## Table of Contents
- [Authorship and Licensing](#authorship-and-licensing)
- [Repository Structure](#repository-structure)
- [Dependencies](#dependencies)
- [Running the Unit Cell Analysis Scripts](#running-the-unit-cell-analysis-scripts)
- [Performing a Non-linear ESL Analysis in Abaqus](#performing-a-non-linear-esl-analysis-in-abaqus)

## Authorship and Licensing

This work includes modified code originally developed by Reinaldo Gonçalves et al. [3]. The original UGENS subroutine script is used with permission from one of the copyright holders. All modifications and extensions are my own.

The remainder of the code in this repository is solely authored by me as a researcher at Aalto University. The GNU General Public License v3.0 applies to all original components. This license permits use, modification, and distribution under the terms outlined therein.

This software is provided "as is", without warranty of any kind. I make no guarantees regarding correctness, performance, or suitability for any particular purpose. You use this code at your own risk. I accept no liability for any errors, data loss, or damages resulting from its use.

## Repository Structure

The files are divided into directories:

- `3D-FEM_Benchmark`    — Abaqus input files used for simulating the benchmark cases.
- `ESL_Analyses`        — Abaqus input files for running the non-linear ESL cases.
- `User_Subroutines`    — Fortran source files defining the calculation of section forces.
- `UC_Analyses`         — Python scripts and additional files for running the unit cell analyses.
- `README.md`           — Project overview and licensing.
- `LICENSE`             — GNU General Public License v3.0.

The Abaqus input files names identify which case the input file is used for in [1]. Example:

```sh
GNL-ESL-R025-2-1-S

GNL-ESL - Geometrically Non Linear Equivalent Single Layer
R025    - Mesh refinement ratio 0.250
2       - Case 2
1       - Load Direction 1 (Compression)
S       - Symmetric Matrix Solver
```

## Dependencies

Running the unit cell analysis requires the following software to be installed:

- Abaqus/CAE.
- Python.

To run the non-linear ESL analyses, you must have:

- Abaqus/CAE.
- A Fortran compiler.
- A software linking the compiler and Abaqus.

These dependencies are proprietary and not included in this repository. You must obtain them separately under their respective licenses.

The code in this repository has been tested with the following software and versions:

- Dassault Systèmes SIMULIA Abaqus/CAE 2023.HF4
- Microsoft Visual Studio Enterprise 2022 v.17.11.5
- Intel Fortran Compiler 2024.1
- Python 3.13

## Running the Unit Cell Analysis Scripts

1. **Open a new Abaqus/CAE environment in the unit cell work directory:** 
The directory must contain

- `parameters.txt`                - specifies the unit cell dimensions, material properties and additional info.
- `UnitCellDistortionMapping.py`  - defines how the ideal unit cell mesh is to be distorted.
- `UnitCellAnalysis.py`           - main python script to be run via Abaqus/CAE.
- `UnitCellPostProcessor.py`      - calculates the stiffness from obtained load response data.

Note that the repository includes two unit cell analysis directories, one for the proposed unit cell analysis and one for the first buckling mode analysis.

2. Input the unit cell parameters in `parameters.txt`:
```sh
- panelName         - string specifying the suffix added to the output files created by "UnitCellPostProcessor.py"
- s  = 0.650000     - unit cell width [mm]
- l  = 0.650000     - unit cell length [mm]
- tp = 0.005000     - plate thickness [mm]
- tw = 0.005000     - stiffener web thickness [mm]
- tf = 0.008870     - stiffener flange thickness [mm]
- hw = 0.091130     - stiffener web height [mm]
- bf = 0.022110     - stiffener flange width [mm]
- ef = 0.011055     - stiffener flange excentricity [mm]
- E  = 206e+09      - panel elastic modulus [Pa]
- nu = 0.3          - panel Poisson's ratio
- nx = 60           - number of elements in x-direction
- ny = 60           - number of elements in y-direction
- M  = 3            - size of distortion definition matrix H
- N  = 3            - size of distortion definition matrix H
```

3. **Define unit cell distortion in `UnitCellDistortionMapping.py`:**
By modifying the H-matrix, you can specify which wavelengths and amplitudes to add to the unit cell model. 

4. **Run `UnitCellAnalysis.py` through Abaqus/CAE by selecting File/Run Script.**
The script will create 12 unit cell models with appropriate boundary conditions and run two jobs (compression and tension of the same load case) at a time. The number of simultaneous jobs is limited to avoid exceeding token limits. The unit cell analysis is completed when the Message "Unit cell analysis completed in X minutes and Y seconds." appears in the Abaqus message window. The generated files `loadCase1.csv` through `loadCase6.csv` contains the extracted data from the 12 analyses, combined into 6 based on the load case.

5. **Run `UnitCellPostProcessor.py` from the command window, PowerShell or similar.**
The script will interpolate the load responses to pre-defined strain points, calculate the stiffness properties, and output 

- `stiffnessMatrix_panelName.f`     - Contains the non-linear stiffness matrix to be copied to the UGENS subroutine.
- `sectionDefinition_panelName.txt`  - Contains the general shell section definition to be copied to the Abaqus Analysis input file.
- `loadResponse.pdf`                 - Plot of load responses as functions of driving strain components.
- `stiffness.pdf`                    - Non-linear stiffness matrix presented in graphical form.

## Performing a Non-linear ESL Analysis in Abaqus

1. **Prepare an Abaqus input file of the ESL geometry:** 
Copy the contents from `sectionDefinition_panelName.txt` to the part or instance definition. Remove the `UNSYMM` keyword in the first line if you wish to run the case with a symmetric matrix solver.

Specify the Riks solver. See the provided input files in this repository for appropriate convergence and time step controls.

2. **Copy the contents of `stiffnessMatrix_panelName.f` to `UGENS-NL-ESL.f`.**
The four versions of the script differ only in the stiffness matrix definition.

3. **Run the analysis from the command window by executing the following command in the directory where the analysis file is located:**
```sh
abaqus job=esl_job_name user=UGENS-NL-ESL.f
```

Note that the relative path to the user subroutine file must be included in the command if the file is located outside of the analysis folder.

To ensure stability of the method, please see the included input files and [2, 3].

## References

[1] Hagnäs, D., Remes, H., Mancini, F and Romanoff, J. (under review). Analysis of Distorted Stiffened Panels under Uniaxial Loading Using Non-linear Equivalent Single Layer Formulation. \

[2] Hagnäs, D. (2025). Modelling of stiffened panels with local distortions using equivalent single layer method [master's thesis]. Espoo, Finland. Aalto University School of Engineering. Available at: https://urn.fi/URN:NBN:fi:aalto-202506094328 \

[3] Reinaldo Gonçalves, B., Jelovica, J. and Romanoff, J. (2016). Abaqus UGENS subroutine for nonlinear analysis of periodic panels. Aalto University publication series SCIENCE + TECHNOLOGY 9/2016. Helsinki, Finland, Aalto University School of Engineering. Available at: https://urn.fi/URN:ISBN:978-952-60-6905-0
