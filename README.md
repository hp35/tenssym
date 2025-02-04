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

## Atlas of non-zero polar and axial tensor elements

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

## Example of output from the `tenssym` program

(From summary: [atlas/rank-3/axial/tensor-hexagonal_6bar-rank-3-axial.txt](atlas/rank-3/axial/tensor-hexagonal_6bar-rank-3-axial.txt)):

```text
Summary generated Tue Feb  4 02:12:53 PM CET 2025 by TensSym for
point-symmetry group hexagonal_6bar, axial tensor of rank 3.
See https://github.com/hp35/tenssym for details.
================================================================================
Symmetry operations associated with 6 (hexagonal, inversion symmetry) are:

1. 6-FOLD ROTATION SYMMETRY AROUND Z-AXIS WITH INVERSION APPLIED ALONG Z AXIS
6z (6-fold rotation symmetry around z-axis with inversion applied along z axis),
described by R(6z) for which det(R(6z)) = -1 (being an improper rotation):
⎡      √3     ⎤
⎢1/2   ──   0 ⎥
⎢      2      ⎥
⎢             ⎥
⎢-√3          ⎥
⎢────  1/2  0 ⎥
⎢ 2           ⎥
⎢             ⎥
⎣ 0     0   -1⎦

Compiling equation system from the 1 symmetries of point-symmetry group 6
(hexagonal, inversion symmetry).
Solving system of 26 equations for nonzero elements.
================================================================================
POINT SYMMETRY GROUP 6 (HEXAGONAL, INVERSION SYMMETRY)

For the axial (pseudo) tensor of rank 3 (for which there in a completely non-
symmetrical case would be a maximum of 27 elements) under constraint of point-
symmetry group 6 (hexagonal, inversion symmetry), there are 21 nonzero elements,
of which 9 are independent.

The 9 independent and nonzero tensor elements are as follows:
       1.  xxx = -xyy = -yxy = -yyx
       2.  xxy = -yyy = xyx = yxx
       3.  xxz = yyz
       4.  xyz = -yxz
       5.  xzx = yzy
       6.  xzy = -yzx
       7.  zxx = zyy
       8.  zxy = -zyx
       9.  zzz = independent
================================================================================
```

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
