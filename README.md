MPIFD-MIMO is a Matlab package that provides a set of tools to identify modal parameters of a structure using measured frequency response functions (FRFs):
- frequency,
- modal damping,
- modal residues,
- real normal mode shapes.

The identification of poles is based on the least-squares complex frequency-domain estimator (LSCF) using fast-stabilizing frequency domain parameter estimation method. The extraction of the stable poles is based on a clear stabilization chart included in an interactive graphical user interface (GUI):

<img width="512" height="418" alt="gui_stabchart" src="https://github.com/user-attachments/assets/2156cfbf-f5cf-418e-b508-2200dc34b87b" />

The residues identification uses the least-square frequency domain estimator (LSFD) and can be performed using classical or non classical damping assumption. Estimated FRFs and mode shapes can be finally viewed in the GUI after extracting and normalizing real normal modes from residues.

The use of the toolbox is detailed step by step in the Matlab notebooks :
- MPIFD_MIMO_4DOF_example.mlx : example of use on a four degrees of freedom numerical system.
- MPIFD_MIMO_plate_example.mlx : example of use on an experimental plate.
