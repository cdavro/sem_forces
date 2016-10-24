## Synopsis

A small fortran (f90) program to calculate force constans (bonds and bend angles force constants for now) in the Amber format using Gaussian09 fchk parsed cartesian Hessian matrix. See below for format.

## Installation

Compile it with gfortran (>= 4.9) or ifort (>= 2013).
Ex: gfortran -framwork Accelerate -o sem_forces.f90 sem_forces (using the MacOS default installed LAPACK framework)

## Usage

Just lauch the program with the data file in the directory (Subject to change and will ask where the data file is located)


## Data file format

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

Here, an example for the H2O molecule (included as dataWater file)

1:  3
2:  --------- Coordinates ---------
3:  1;8;8.62816615E-32;-3.60822483E-16;2.09404098E-01
4:  2;1;-1.97215226E-31;1.48149974E+00;-8.37616393E-01
5:  3;1;-1.81431391E-16;-1.48149974E+00;-8.37616393E-01
6:  --------- Constants ---------
7:  1;1;1;2;2;2;3;3;3
8:  8;8;8;1;1;1;1;1;1
9:  O;O;O;H;H;H;H;H;H
10:  x1;y1;z1;x2;y2;z2;x3;y3;z3
11: -1.64542770E-02;-2.54708991E-13;1.60142221E-12;8.22713851E-03;5.47966221E-13;8.85821784E-13;8.22713851E-03;-1.04864019E-14;-1.69199007E-12
12: -2.54708991E-13;7.43613112E-01;-2.85215185E-11;-1.13640799E-12;-3.71806556E-01;2.68616207E-01;7.68294756E-13;-3.71806556E-01;-2.68616207E-01
13: 1.60142221E-12;-2.85215185E-11;3.84979958E-01;-2.17658745E-12;1.93033237E-01;-1.92489979E-01;1.38171120E-12;-1.93033237E-01;-1.92489979E-01
14: 8.22713851E-03;-1.13640799E-12;-2.17658745E-12;-7.99167272E-03;-9.44132586E-14;5.77755922E-13;-2.35465793E-04;1.37388978E-13;2.47355305E-12
15: 5.47966221E-13;-3.71806556E-01;1.93033237E-01;-9.44132586E-14;4.02401599E-01;-2.30824722E-01;1.37386376E-13;-3.05950427E-02;3.77914855E-02
16: 8.85821784E-13;2.68616207E-01;-1.92489979E-01;5.77755922E-13;-2.30824722E-01;1.87309946E-01;-2.47354961E-12;-3.77914855E-02;5.18003321E-03
17: 8.22713851E-03;7.68294756E-13;1.38171120E-12;-2.35465793E-04;1.37386376E-13;-2.47354961E-12;-7.99167272E-03;-9.50196847E-14;-5.78137039E-13
18: -1.04864019E-14;-3.71806556E-01;-1.93033237E-01;1.37388978E-13;-3.05950427E-02;-3.77914855E-02;-9.50196847E-14;4.02401599E-01;2.30824722E-01
19: -1.69199007E-12;-2.68616207E-01;-1.92489979E-01;2.47355305E-12;3.77914855E-02;5.18003321E-03;-5.78137039E-13;2.30824722E-01;1.87309946E-01


You can write a fchk parser
The parsed file should be called 'data'

## License

Distributed under the MIT license
