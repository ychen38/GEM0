# GEM - Software for kinetic simulation of magnetically confined fusion plasmas
The Direct directory contains source files for the split-weight control-variate mathod

The Hybrid directory contains source files for the gyrokinetic ion/fluid electron hybrid model

Good Software Engineering Guidelines:

1, Break long lines into short continued lines (with &)
2, No tabs. Use whitespaces.  
3, All reals should be promoted by compiler; so, 'REAL(8)' in the code should just be 'real'
4, No old-style numbered "do/continue" loops. Use "do/end do" with no numbers or "continue"
Old-style "character*70" should be changed to "character(len=70)"
5, "Return" is typically not needed in subroutines.
6, Subroutines should end with 'end subroutine myroutine', not just 'end'.
7, Try to convert all .f files to .f90 with free-format syntax
8, All declared variables should be used. Remove variables not needed anymore. 
