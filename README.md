# Computation of nonzero elements of polar and axial tensors by direct inspection

Computation of nonzero and independent elements of arbitrary-rank polar or
axial tensors ("true tensors" or "pseudotensors", respectively). This is done
symbolically by applying a scheme of direct inspection for the symmetry
operations allowed by the targeted point-symmetry group, and is handled by
the symbolic math library SymPy of Python. See https://docs.sympy.org/

For the analysis of symmetry within the crystallographic point-symmetry groups,
the very instructive tutorial available at the University of Oklahoma, Dept of
Chemistry and Biochemistry, at http://xrayweb.chem.ou.edu/notes/symmetry.html
can be recommended as a starter.

References

  1. Susceptibility Tensors for Nonlinear Optics, S. V. Popov, Yu. P. Svirko
     and N. Zheludev (Institute of Physics Publishing, Bristol, 1995).
     ISBN 0-7503-0253-4.

  2. The Elements of Nonlinear Optics, P. N. Butcher and D. Cotter
     (Cambridge University Press, Bristol, 1993). ISBN 0-521-42424-0.

  3. The Nonlinear Optics of Magneto-Optic Media (doctoral thesis),
     Fredrik Jonsson (The Royal Institute of Technology, 2000),
     ISBN 91-7170-575-9; Available via the National Library of Sweden (Kungl.
     Biblioteket) at http://urn.kb.se/resolve?urn=urn:nbn:se:kth:diva-2967

## Compiling and running the code

The `tenssym.sh` script has been used to generate the enclosed [atlas](./atlas/)
of nonzero tensor elements for the following point-symmetry groups (tensor
ranks from two to five):

    "cubic_4bar3m"   - 43m cubic, inversion symmetry
    "cubic_23"       - 23 cubic, no inversion symmetry
    "trigonal_3"     - 3 trigonal, no inversion symmetry
    "trigonal_32_2x" - 32 trigonal, no inversion symmetry, 2-fold x-symmetry
    "trigonal_32_2y" - 32 trigonal, no inversion symmetry, 2-fold y-symmetry
    "triclinic_1bar" - 1 triclinic, inversion symmetry for all axes
    "triclinic_1"    - 1 triclinic, no inversion symmetry
    "tetragonal_422" - 422 tetragonal, no inversion symmetry
    "cubic_m3m"      - m3m cubic, inversion symmetry
    "cubic_432"      - 432 cubic, no inversion symmetry
    "cubic_m3"       - m3 cubic, inversion symmetry
    "hexagonal_6"    - 6 hexagonal, no inversion symmetry
    "hexagonal_6bar" - 6 hexagonal, inversion symmetry

## Compiling and running the code

The code can be run from the enclosed script `tenssym.sh`, iterating over a
set of implemented point-symmetries and generating the entire catalogue of
non-zero tensor elements.

```bash
$ cd tenssym
$ ./tenssym.sh
```

## Copyright
Copyright (C) 2023, Fredrik Jonsson, under GPL 3.0. See enclosed LICENSE.

## Location of master source code
The source and documentation can be found at https://github.com/hp35/tenssym
