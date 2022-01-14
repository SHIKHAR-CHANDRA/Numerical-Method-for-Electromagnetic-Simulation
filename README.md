# Numerical-Method-for-Electromagnetic-Simulation

A 1-D Finite Difference Time Domain algorithm is used to simulate electromagnetic wave behavior.

The governing equations for the algorithm to work are the Maxwell's equations.

The curl equations of Electric and Magnetic Fields are staggered both in time for finite difference approximation.

The Yee Grid scheme allows the electric and magnetic fields' spatial components to be staggered in space in order to have a divergence free grid.

Update equations for both the magnetic and electric field are run over a loop in a 'Leap-Frog Algorithm'. This loop runs over all values of position in the grid.

Using Dirichlet Boundary conditions and assuming the fields outside the grid are 0, the fields are updated at the start and end positions.

Grid resolution, or, the step size in space, is calculated using the minimum wavelength in the grid. 

Using Courant Stability Condition, which states that the physical wave cannot propagate farther than a single unit cell in one time step, time resolution is calculated.

The source is a Gaussian pulse which is introduced after small delay. 

The grid is partitioned into Total Field Region and Scattered Field Region.



