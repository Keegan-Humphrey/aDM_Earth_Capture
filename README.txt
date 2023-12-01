Author: Keegan Humphrey
Creation Date: 28/07/23

This repository was created to calculate the Energy Loss, Capture, and Equilibrium Configuration of atomic Dark Matter (aDM) in the Earth over a wide range of masses, velocities, kinetic mixings, and abundances. 

Function definitions are organized into packages in the /package/ folder. One can either open each and click "run all code" according to the "Needs" commands of the .nb file being used, or permanently load them into Mathematica. 

The pipeline to determine the Energy Loss and minimum kinetic mixing needed for capture, is to 
1) run the notebook in /extract optical data/, which fits a sum of Drude oscillators to measured optical data, extracting and saving parameters which define a data driven Mermin Dielectric for the material
2) run the notebook in /compute energy loss/, which computes a grid of energy loss (in [eV m^-1]) over a range of kinematics, and determines an interpolating function as a function of the kinematics saving them to an association. This is the most computationally intensive part of the pipeline
3) run the notebook in /kappa min/, which computes the minimum kinetic mixing need for capture as a function of mass and velocity for a simple model of the Earth using the computed energy loss interpolation

Currently, this is based exclusively on interaction with SM electrons. We look to include nuclear interactions shortly. 

We also look to determine the equilibrium configuration using a simple Vlasov Poisson model of the captured DM, which we intend to solve numerically.