# bonCar
A Fortran code of converting bond lengths to Cartersian coordinates. 

# Introduction
Converting Cartesian coordinates to bond lengths is simpler than the inverse convert. The bond-Cart transformation is intuitively straight-forward but hard to realize. 

## Difficulty
The first difficulty I met is how to solve the fourth coordinates if we knew three coordinates and three lengths, which confused me for many days. Actually this is a common transformation in GIS called [**trilateration**](https://math.stackexchange.com/questions/2969363/finding-a-4th-point-in-3d-space-knowing-3-other-points-and-2-distances-to-the-4t). 

The second difficulty I didn't solve is how to deal with the colinear condition. If three or more atoms are in one line, there will not be a solution for the fourth point, although the condition is rare. 

## Converting Procedure
1. The first atom was placed at the origin. 
2. The second atom was placed at the x axis. 
3. The third atom was placed at the xOy plane. 
4. The fourth was placed at the positive half Z axis. 
5. Then other atoms was solved one by one. 

# Usage
The core subroutine is `bonCar`. The number of bond lengths, bond lengths, the number of atoms are required. 

Intel Fortran Compiler is needed. No more prerequisite. 

1. Clone this repository. 
2. `make`
