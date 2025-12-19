## Overview

MPIFD-MIMO is a Matlab package that provides a set of tools to identify modal parameters of a structure using measured Frequency Response Functions (FRFs):
- frequency,
- modal damping,
- modal residues,
- real normal mode shapes.

The identification of poles is based on the least-squares complex frequency-domain estimator (LSCF) using fast-stabilizing frequency domain parameter estimation method. The extraction of the stable poles is based on a clear stabilization chart included in an interactive graphical user interface (GUI). The residues identification uses the least-square frequency domain estimator (LSFD) and can be performed using classical or non classical damping assumption. Estimated FRFs and mode shapes can be finally viewed in the GUI after extracting and normalizing real normal modes from residues. An overview of the identification process can be described as follow:

<img width="661.5" height="618" alt="figure" src="https://github.com/user-attachments/assets/ffacb588-edf0-45e8-b9b6-a25a0d518c90" />

## Example files

The use of the toolbox is detailed step by step in the Matlab notebooks:
- _MPIFD_MIMO_4DOF_example.mlx_ : example of use on a four degrees of freedom numerical system.
- _MPIFD_MIMO_plate_example.mlx_ : example of use on an experimental plate.

Associated examples can be directly runned without notebook:
- _4dof_example/four_dof_example_script.m_
- _plate_example/plate_example_script.m_

## Reference

Chomette et al. (2025). MPIFD-MIMO : A Matlab tool for modal parameters identification in frequency domain. _Journal of Open Source Software_, VOL (ISSUE), PAGE, https://doi.org/10.xxxxxx/draft.
