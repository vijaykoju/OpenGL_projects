------------------Andromeda-MilkyWay Collision Simulation--------------------

Programmer  : Vijay Koju

-----------------------------------------------------------------------------

This program simulates collision dynamics of Andromeda and MilkyWay galaxies. This kind of problem is often called N-Body problem, and in general the complexity of such a problem is O(n^2) for each simulation step because it requires the computation of gravitational force on each body by all other bodies. There are many other algorithm such as Barnes-Huts Tree algorithm (O(nlogn)), and Fast Multipole Method (O(n)) which use some approximations to reduce the time complexity. This program, however, uses none of such algorithms. It approximates the two galaxies as two stars located at the center of mass of the galaxies with some huge mass. These two stars interact according to the laws of gravitation. All the stars in both the galaxies are only influenced by these two core stars. The numerical integration in this program is done using the fourth order Hermite Scheme.

-----------------------------------------------------------------------------

******************* Required files to run this program **********************
1) FinalPro_GalSim.cpp
2) intialConditionData.dat
3) CoolWarmFloat257.dat
4) Jet257.dat
5) makefile (for a MAC OS)

****** To run this program (if all the files listed above are present) ******
    $ make
    $ ./galSim

--------------------Keyboard options:----------------------

h                   --> view this page                ##

s                   --> start/stop animation          ##

f                   --> fullscreen mode on/off        ##

c                   --> change color warm/cool        ##

o                   --> boundary (box) on/off         ##

i                   --> go back to initial setting    ##

r                   --> reset with different initial velocity                      ##

t                   --> trace the path of "Center of mass" of the two galaxies     ##

q or esc            --> quit the program              ##

Ctrl +              --> zoom in                       ##

- (minus sign)      --> zoom out                      ##
- 
 Right arrow         --> move to the right             ##

Left arrow          --> move to the left              ##

Up arrow            --> move up                       ##

Down arrow          --> move down                     ##

-----------------------Mouse options:----------------------

Left click and drag --> rotate                        ##

Right click         --> pop menu                      ##

-----------------------------------------------------------------------------
