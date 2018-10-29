# sfoglie
Coupled Boundary Layer Theory / Potential Theory Boundary Element Method solver   
for the incompressible flow around an airfoil

![Profile Image](https://github.com/davecats/sfoglie/blob/master/profile_img.png)

*sfoglie* is a MATLAB port of the subsonic airfoil development system [XFOIL](https://web.mit.edu/drela/Public/web/xfoil/) 
by M. Drela, originally written in F77. *sfoglie* ports most (yet not all) of the XFOIL functionalities and adds a few, most
notably the possibility to add a wall-normal velocity distribution along the profile chord, which simulates flow control
by blowing or suction. One missing features compared to XFOIL is the solution of the aerodynamic inverse problem.


