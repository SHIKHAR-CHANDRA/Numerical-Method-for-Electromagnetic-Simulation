# Numerical-Method-for-Electromagnetic-Simulation

A 1-D Finite Difference Time Domain algorithm is used to simulate electromagnetic wave behavior.

The governing equations for the algorithm to work are the Maxwell's equations.

The curl equations of Electric and Magnetic Fields are staggered in time for finite difference approximation.

The Yee Grid scheme allows the electric and magnetic fields' spatial components to be staggered in space in order to have a divergence free grid.

Update equations for both the magnetic and electric field are run over a loop in a 'Leap-Frog Algorithm'. This loop runs over all values of position in the grid.

Using Dirichlet Boundary conditions and assuming the fields outside the grid are 0, the fields are updated at the start and end positions.

Grid resolution, or, the step size in space, is calculated using the minimum wavelength in the grid. 

Using Courant Stability Condition, which states that the physical wave cannot propagate farther than a single unit cell in one time step, the time resolution is calculated.

The source is a Gaussian pulse and is injected as a soft source after a small delay.

The grid is partitioned into Total Field Region and Scattered Field Region at the point where the source is injected. The scattered field region contains only the backward scattered wave field components whereas the total field region contains the transmitted as well as scattered waves from the materials. 

The complete algorithm is run over a loop in time and each iteration at a time step solves the update equations at all points in the grid.

To calculate the results, i.e, the reflectance and transmittance, Discrete Time Fourier Transform is performed. 


# Radome Design

The problem is to design a radome to protect an antenna of certain width and operating at specified frequency. The dielectric constant of the material of antenna is given and the goal is to maximize transmission through the radome (have negligible reflectance). 

For this problem, anti-reflection layers are added to both sides of the antenna in the form of quarter wave transformer structures. The dielectric constant and refractive index of the radome is calculated, using which, the thickness of the radome is determined. 

The materials, namely, the antenna and two radome layers, are added in the 1-D grid by changing the values of relative permittivity at the respective grid positions. 

Simulating the FDTD algorithm, the results show that the reflectance is negligible at the operating frequency.



