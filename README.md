# Fix spring/self/move

## Description
This fix creates a virtual position for every atom it is invoked upon,
and tethers the atom to the virtual position by a 1-D spring. The virtual
position is moved by the amount specified in the fix every timestep. The atom position is
also moved by the amount specifed in the fix, but this move is in addition 
to the position update due to other fixes. 

<img src="fix_spring_working.svg" width="400" height="200" />

The fix accepts a numeric value for the spring constant, and [equal style
variables](https://docs.lammps.org/variable.html) (or NULL) for the position update values. 

The virtual position is stored in the variable 'xoriginal'. It is
initialized to be equal to the atom position when the fix is invoked.

## Compilation: 
To use this fix, place fix_spring_move.cpp and fix_spring_move.h 
in the LAMMPS src directory, and compile as usual. This fix has been tested
with the 29 SEPT 2021 version of LAMMPS.

## Usage: 
```
fix springm walls   spring/self/move 10                 xyz                v_smov NULL NULL
```
fix fixID   groupID fix_name         spring_constant    move_directions    equal_style_variables

This fix can be used to allow a group of atoms maintain its general shape
while being moved, with each individual atom also allowed to have thermal
motion.
