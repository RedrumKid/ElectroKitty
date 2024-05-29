ElectroKitty is a simulator built for simulating different mechanisms on an RDE electrode geometry

Repository containing the main simulator code and some examples on how to run it
Examples are given in the Tutorial folder to get a grasp on how the simulator runs

The program runs by running the ElectroKitty_main.py script, so have that in the same folder as the script running it

The simulation requieres that you declare a mechanism as a string (check the tutorial for details),
and based on these, given the parameters, which the user provides as lists of constants calculates the appropriate current response,
given the potential program. 

Currently the simulator can handle these techniques: CV, ACV, and CA. Other potentiodynamic techniques can be implemented,
with the user supplying a pregenerated signal, but they have not been tested yet.

Version history:
(currently the simulator is in an alpha, so some functionalites may not work properly)

Alpha 1.1:  Now with added roatation and code cleanup

Alpha 1.2: Kinetic list made intuitive

Alpha 1.3: Implemented multi-elecron support for BV kinetics

Alpha 1.4: Implementecd isotherm modeling