# How to run the membrane code

The code needs to be run in a folder containing an input script with parameters, 'in.ves' as well as an initial configuration scripts: eg 'icos_5882_3_5.dat' (the numbers are: number of membrane beads, number of colloids and the diameter of the colloids)

## Python scripts for making initial configuration and the input script for the c-code

Some example python scripts are contained in the folder 'Python_MakingInputs' which takes parameters eg:

build_tether.py ${arc1} ${arc2} ${tau} ${sigma} ${D0} ${seed}

where arc1 and arc2 are the two arc distances between the 3 colloids and tau is the angle that these two arcs make when they meet. Sigma is the diameters of the colloids, D0 is the interaction energy and the seed dictates the random initial velocity.

## Compile the code:

Compile code with:
gcc Vescicle_Tether.c -o Vescicle_Tether -lm

to the run with:
./Vescicle_Tether in.ves

if in the same folder as in initial data file with the name icos_5882_3_5.dat