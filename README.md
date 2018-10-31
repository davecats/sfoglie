# sfoglie
Coupled Boundary Layer Theory / Potential Theory Boundary Element Method solver   
for the incompressible flow around an airfoil

![Profile Image](https://github.com/davecats/sfoglie/blob/master/profile_img.png)

*sfoglie* is a MATLAB port of the subsonic airfoil development system [XFOIL](https://web.mit.edu/drela/Public/web/xfoil/) 
by M. Drela, originally written in F77. *sfoglie* ports most (yet not all) of the XFOIL functionalities and adds a few, most
notably the possibility to add a wall-normal velocity distribution along the profile chord, which simulates flow control
by blowing or suction. One missing features compared to XFOIL is the solution of the aerodynamic inverse problem.

### Usage

The file documentation.pdf contains a description of the program and usage instruction. If you are in a hurry and want to start quickly, just run *main.m* in MATLAB and see the output!  

### Contacts

Dr. Davide Gatti  
davide.gatti [at] kit.edu  
msc.davide.gatti [at] gmail.com  

Karlsruhe Institute of Technology  
Institute of Fluid Dynamics  
Kaiserstraße 10  
76131 Karlsruhe  

### How to cite this code

If you use this code and find it helpful, please cite:  
``` M.Reder, A. Stroh and D. Gatti, "Preliminary study of flow control via uniform blowing on airfoils with a boundary element method", Notes on Numerical Fluid Mechanics and Multidisciplinary Design, (submitted, 2018) ```
