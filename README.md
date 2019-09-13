# Virtual_Material_Lab
The code simulates and visualizes displacement controlled mechanical tests on a unit sized material. Single scale and multiscale analyses can be performed. The visualizations are then saved as gifs.   

The virtual material lab essentially allows the user to see the effects of material properties on purely mechanical tests. The multiscale analyses allows the visualization of two phased composites (mostly fiber and matrix) and additionally shows the ffects of matrix and fiber degradation.

The code solves a linear displacement field such that it satisfies the equilibrium, continuity, and constiutive equations. The solutions to these equations leads to values for the displacement field parameters. The displacement field is discretized into 9 nodes over a square (unit material) and the displacement is visualized on application of a strain field. The strain field can be uniaxial, biaxial or triaxial and a combination of shear and normal strains, although care has to be taken to avoid overconstraining the system. This is performed by a smart reordering of the constitutive equations so that all the knowns and unknowns are consistent. See the choice option of [0 0 0 1] which implements this.

Numerical differentiation is performed using the forward euler method.   
