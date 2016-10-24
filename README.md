## What is it ?

A small fortran (f90) program to calculate force constans (bonds and bend angles force constants for now) in the Amber format using Gaussian09 fchk parsed cartesian Hessian matrix. See below for format.  

## Installation

Compile it with gfortran (>= 4.9) or ifort (>= 2013).  
Ex: gfortran -framwork Accelerate -o sem_forces.f90 sem_forces (using the MacOS default installed LAPACK framework)  

## Usage

Just lauch the program with the data file in the directory (Subject to change and will ask where the data file is located)  


## Data file format

You can write a fchk parser.  
The parsed file should be called 'data'.  

1: Number of Atoms  
2: --------- Coordinates --------- (or blank line)  
3: Number1;AtomicNumber1;X1;Y1;Z1  (in bohr)
4: Number2;AtomicNumber2;X2;Y2;Z2
...
.: --------- Constants --------- (or blank line)
.: Number1;Number1;Number1;Number2;Number2;Number2;...  
.: AtomicNumber1;AtomicNumber1;AtomicNumber1;AtomicNumber2;AtomicNumber2;AtomicNumber2;...  
.: AtomicSymbol1;AtomicSymbol1;AtomicSymbol1;AtomicSymbol2;AtomicSymbol2;AtomicSymbol2;...  
.: x1;y1;z1;x2;y2;z2;...  
.: Fx1x1;Fy1x1;Fz1x1;Fx2x1;...   (in hartree/bohr^2)  
.: Fx1y1;Fy1y1;Fz1y1;Fx2y1;...  
.: Fx1z1;Fy1z1;Fz1z1;Fx2z1;...  
...  

An example for the H2O molecule is included (as dataWater file, rename it to 'data' to use it)  

## License

Distributed under the MIT license  
