This repository was created to calculate the Energy Loss, Capture, and Equilibrium Configuration of atomic Dark Matter (aDM) in the Earth over a wide range of masses, velocities, kinetic mixings, and abundances. The documentation folder contains a set of detailed notes with references that were used to construct this pipeline.

Function definitions are organized into packages in the /package/ folder. One can either open each and click "run all code" according to the "Needs" commands of the .nb file being used, or permanently load them into Mathematica. 

The pipeline to determine the Energy Loss and minimum kinetic mixing needed for capture, is to 
1) run the notebook in /extract optical data/, which fits a sum of Drude oscillators to measured optical data, extracting and saving parameters which define a data driven Mermin Dielectric for the material
2) run the notebook in /compute energy loss/, which computes a grid of energy loss (in [eV m^-1]) over a range of kinematics, and determines an interpolating function as a function of the kinematics saving them to an association. This is the most computationally intensive part of the pipeline
3) run the notebook in /kappa min/, which computes the minimum kinetic mixing need for capture as a function of mass and velocity for a simple model of the Earth using the computed energy loss interpolation

Modules to compute nuclear scattering are defined in EnergyLoss.nb. There are both 0 temperature and finite temperature modules (both in the non-degenerate / classical limit). 

The capture folder contains a notebook to implement a capture algorithm described in the notes. At intermediate interaction strengths, this consists of a detailed monte-carlo propagation algorithm. 

We also look to determine the equilibrium configuration using a simple Vlasov Poisson model of the captured DM, which we intend to solve numerically.