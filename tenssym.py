#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
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

 [1] Susceptibility Tensors for Nonlinear Optics, S. V. Popov, Yu. P. Svirko
     and N. Zheludev (Institute of Physics Publishing, Bristol, 1995).
     ISBN 0-7503-0253-4.

 [2] The Elements of Nonlinear Optics, P. N. Butcher and D. Cotter
     (Cambridge University Press, Bristol, 1993). ISBN 0-521-42424-0.

 [3] The Nonlinear Optics of Magneto-Optic Media (doctoral thesis),
     Fredrik Jonsson (The Royal Institute of Technology, 2000),
     ISBN 91-7170-575-9; Available via the National Library of Sweden (Kungl.
     Biblioteket) at http://urn.kb.se/resolve?urn=urn:nbn:se:kth:diva-2967

Copyright (C) 2023, Fredrik Jonsson.
"""
import math, sys, getopt
import numpy as np
import sympy as sp
import textwrap
from enum import Enum
from typing import List, Tuple
from dataclasses import dataclass
from sympy import init_printing
from collections import Counter

init_printing()

def shannon(string):
    """
    Compute the Shannon entropy of a string of characters. For a string of all
    characters being equal (even for an empty string), the Shannon entropy is
    zero, with higher entropies indicating a higher "mixing ratio".
    
    When computing non-zero tensor elements from symmetries, we in automatic
    mode use the Shannon entropy to sort out the elements which preferrably
    should be chosen as the unique ones. (Typically, for a rank-4 tensor,
    this set would start with "xxxx", "yyyy", "zzzz", followed by "xxyy",
    "xyxy" etc.)

    Examples of entropies:
        "xx", "yy", "xxxx", "yyyy", "zzzz":   Entropy = 0.0
        "xy", "yx", "xxyy", "xyyx", "xyxy":   Entropy = 1.0
        "xyyy", "yyyx", "zyyy", "yyyz":       Entropy = 0.81
        "xyz", "xyzxyz", "xyzxyzxyz" :        Entropy = 1.58

    See https://en.wikipedia.org/wiki/Entropy_(information_theory)

    Parameters
    ----------
    string : str
        The string of characters which we wish to compute the Shannon
        entropy for.

    Returns
    -------
    float
        The Shannon entropy of the supplied string.
    """
    counts = Counter(string)
    frequencies = ((i/len(string)) for i in counts.values())
    return -sum(f*math.log2(f) for f in frequencies)

""" Tensor type (polar or axial) """
class TensorType(Enum):
    POLAR = "polar (true)"   # "True tensor"
    AXIAL = "axial (pseudo)" # "Pseudotensor"

""" Crystallographic point-symmetry groups """
class PointSymmetryGroup(Enum):
    cubic_4bar3m   = "43m (cubic, inversion symmetry)"
    cubic_23       = "23 (cubic, no inversion symmetry)"
    cubic_m3m      = "m3m (cubic, inversion symmetry)"
    cubic_432      = "432 (cubic, no inversion symmetry)"
    cubic_m3       = "m3 (cubic, inversion symmetry)"
    trigonal_3     = "3 (trigonal, no inversion symmetry)"
    trigonal_32_2x = "32 (trigonal, no inversion symmetry, 2-fold x-symmetry)"
    trigonal_32_2y = "32 (trigonal, no inversion symmetry, 2-fold y-symmetry)"
    triclinic_1    = "1 (triclinic, no inversion symmetry)"
    triclinic_1bar = "1 (triclinic, inversion symmetry for all axes)"
    tetragonal_422 = "422 (tetragonal, no inversion symmetry)"
    hexagonal_6    = "6 (hexagonal, no inversion symmetry)"
    hexagonal_6bar = "6 (hexagonal, inversion symmetry)"

""" Axis of rotation for rotation matrices """
class Axis(Enum):
    X = "x"                  # The cartesian x-axis
    Y = "y"                  # The cartesian y-axis
    Z = "z"                  # The cartesian z-axis
    XYZ = "xyz"              # The diagonal (x,y,z)-axis
    NONE = "none"            # No axis

""" Angle of rotation for rotation matrices """
class Angle(Enum):
    ZERO = "0"               # Zero radians   (no rotational symmetry)
    PI = "pi"                # PI radians     (2-fold rotational symmetry)
    TWO_PI_THIRD = "2pi/3"   # 2*PI/3 radians (3-fold rotational symmetry)
    PI_HALF = "pi/2"         # PI/2 radians   (4-fold rotational symmetry)
    PI_THIRD = "pi/3"        # PI/3 radians   (6-fold rotational symmetry)

class RotationMatrix:
    def __init__(self, axis, angle, inv=None):
        self.axis = axis
        self.angle = angle

        if (angle == Angle.ZERO) or (axis == Axis.NONE):
            self.oplabel = "I"
            self.opdescr = "identity operation"
            self.array = sp.Matrix([
                [  1,  0,  0 ],
                [  0,  1,  0 ],
                [  0,  0,  1 ]])
            if (inv!=None): self.invert(inv)
            return

        if axis == Axis.X:   # Rotation around the x-axis

            if angle == Angle.PI_HALF:
                self.oplabel = "4x"
                self.opdescr = "4-fold rotation symmetry around x-axis"
                s32 = sp.sqrt(3)/2
                r12 = sp.Rational(1,2)
                self.array = sp.Matrix([
                    [  1,  0,  0 ],
                    [  0,  0,  1 ],
                    [  0, -1,  0 ]])
                if (inv!=None): self.invert(inv)
                return

            if angle == Angle.PI:
                self.oplabel = "2x"
                self.opdescr = "2-fold rotation symmetry around x-axis"
                self.array = sp.Matrix([
                    [  1,  0,  0 ],
                    [  0, -1,  0 ],
                    [  0,  0, -1 ]])
                if (inv!=None): self.invert(inv)
                return

        if axis == Axis.Y:   # Rotation around the y-axis

            if angle == Angle.PI_HALF:
                self.oplabel = "4y"
                self.opdescr = "4-fold rotation symmetry around y-axis"
                self.array = sp.Matrix([
                    [  0,  0, -1 ],
                    [  0,  1,  0 ],
                    [ -1,  0,  0 ]])
                if (inv!=None): self.invert(inv)
                return

            if angle == Angle.PI:
                self.oplabel = "2y"
                self.opdescr = "2-fold rotation symmetry around y-axis"
                self.array = sp.Matrix([
                    [ -1,  0,  0 ],
                    [  0,  1,  0 ],
                    [  0,  0, -1 ]])
                if (inv!=None): self.invert(inv)
                return

        if axis == Axis.Z:   # Rotation around the z-axis

            if angle == Angle.PI_THIRD:
                self.oplabel = "6z"
                self.opdescr = "6-fold rotation symmetry around z-axis"
                r12 = sp.Rational(1,2)
                s32 = sp.sqrt(3)/2
                self.array = sp.Matrix([
                    [  r12,  s32,  0 ],
                    [ -s32,  r12,  0 ],
                    [    0,    0,  1 ]])
                if (inv!=None): self.invert(inv)
                return

            if angle == Angle.PI_HALF:
                self.oplabel = "4z"
                self.opdescr = "4-fold rotation symmetry around z-axis"
                self.array = sp.Matrix([
                    [  0,  1,  0 ],
                    [ -1,  0,  0 ],
                    [  0,  0,  1 ]])
                if (inv!=None): self.invert(inv)
                return

            if angle == Angle.TWO_PI_THIRD:
                self.oplabel = "3z"
                self.opdescr = "3-fold rotation symmetry around z-axis"
                r12 = sp.Rational(1,2)
                s32 = sp.sqrt(3)/2
                self.array = sp.Matrix([
                    [ -r12,  s32,  0 ],
                    [ -s32, -r12,  0 ],
                    [    0,    0,  1 ]])
                if (inv!=None): self.invert(inv)
                return

            if angle == Angle.PI:
                self.oplabel = "2z"
                self.opdescr = "2-fold rotation symmetry around z-axis"
                self.array = sp.Matrix([
                    [ -1,  0,  0 ],
                    [  0, -1,  0 ],
                    [  0,  0,  1 ]])
                if (inv!=None): self.invert(inv)
                return

        if axis == Axis.XYZ: # Rotation around the diagonal xyz-axis

            if angle == Angle.TWO_PI_THIRD:
                self.oplabel = "3xyz"
                self.opdescr = "3-fold rotation symmetry around xyz-diagonal"
                self.array = sp.Matrix([
                    [  0,  1,  0 ],
                    [  0,  0,  1 ],
                    [  1,  0,  0 ]])
                if (inv!=None): self.invert(inv)
                return

        raise ValueError("Unimplemented (and unneccessary) rotation %s."%axis)
        return

    def invert(self, invAxis):
        """
        Invert the currently held rotation matrix with respect to the
        supplied axis.

        Parameters
        ----------
        invAxis : Axis
            The axis along which the rotation matrix will be inverted.

        Returns
        -------
        None. The RotationMatrix object will instead be modified accordingly.
        """
        if invAxis == Axis.X: cols = [0]
        elif invAxis == Axis.Y: cols = [1]
        elif invAxis == Axis.Z: cols = [2]
        elif invAxis == Axis.XYZ: cols = [0,1,2]
        else: raise ValueError("Unidentified axis of inversion '%s'."%invAxis)
        for row in range(3):
            for col in cols:
                self.array[row,col] = -self.array[row,col]
        self.opdescr += " with inversion applied along %s axis"%invAxis.value
        return

    def det(self):
        """
        Compute the determinant of the currently contained rotation matrix.
        The computation is done by symbolic computation and yields exact
        results without any floating-point round-off errors.

        Returns
        -------
        float
            The determinant of the rotation matrix. Should always be
            either +1 or -1 (exactly).
        """
        return self.array.det()

    def describesProperRotation(self, tol=1.0e-10):
        """
        Check if the currently contained rotation matrix describes a proper
        (det(R) = 1) or improper (det(R) = -1) rotation (the latter being a
        proper rotation followed by an inversion). Meanwhile, also check that
        the absolute value of the determinant is unity within the supplied
        tolerance; if not, then the contained matrix neither describes a
        proper nor improper rotation, and there is then something seriously
        wrong (for which case an error will be thrown).

        Parameters
        ----------
        tol : float, optional
            Tolerance for the maximum deviation of the absolute value of the
            determinant from unity. Default tolerance is 1.0e-10.

        Returns
        -------
        bool
            True if found to be a proper rotation, for det(R) = 1; otherwise
            False is returned (for an improper rotation including an inversion,
            for which det(R) = -1).
        """
        det = self.det()
        dev = abs(det)-1.0
        if dev > tol:
            raise ValueError("Absolute value of determinant outside of "
                             "tolerance (%e). This should never happen."%dev)
        return (det > 0.0)

    def summarize(self):
        print(textwrap.fill("%s (%s), described by R(%s) for which det(R(%s)) "
              "= %s (being a%s rotation):"%(self.oplabel, self.opdescr,
              self.oplabel,self.oplabel, self.det(),
              " proper" if self.describesProperRotation() else "n improper"),80))
        sp.pprint(self.array)
        return

class PointSymmetryOperations:
    """
    The PointSymmetryOperations class provides tabulated symmetry operations
    for each point-symmetry group.

    For details on the symmetry operations of each group, see Ref.[1],
    Tables 3.3 and 3.5 (pages 41 and 47, respectively).
    """
    def __init__(self, symmetrygroup=PointSymmetryGroup.cubic_23):
        self.symmetrygroup = symmetrygroup

        if (symmetrygroup == PointSymmetryGroup.cubic_23):
            """
            Point-symmetry group 23 (cubic, no inversion symmetry), providing
            symmetry under the three transformations:
                2z     (2-fold rotation symmetry around z, no inversion).
                3xyz   (3-fold rotation symmetry around xyz-diagonal)
                2x     (2-fold rotation symmetry around x, no inversion).
            """
            self.numSymmetryOperations = 3
            self.symmetryOperations = np.array([
                RotationMatrix(Axis.Z, Angle.PI),             # "2z" symmetry
                RotationMatrix(Axis.XYZ, Angle.TWO_PI_THIRD), # "3xyz" symmetry
                RotationMatrix(Axis.X, Angle.PI)])            # "2x" symmetry

        elif (symmetrygroup == PointSymmetryGroup.tetragonal_422):
            """
            Point-symmetry group 422 (tetragonal, no inversion symmetry),
            providing symmetry under the four transformations:
                4z     (4-fold rotation symmetry around z, no inversion)
                2x     (2-fold rotation symmetry around x, no inversion).
                2y     (2-fold rotation symmetry around y, no inversion).
                2z     (2-fold rotation symmetry around z, no inversion).
            """
            self.numSymmetryOperations = 4
            self.symmetryOperations = np.array([
                RotationMatrix(Axis.Z, Angle.PI_HALF), # "4z"
                RotationMatrix(Axis.X, Angle.PI),      # "2x"
                RotationMatrix(Axis.Y, Angle.PI),      # "2y"
                RotationMatrix(Axis.Z, Angle.PI)])     # "2z"

        elif (symmetrygroup == PointSymmetryGroup.cubic_4bar3m):
            """
            Point-symmetry group 43m (cubic, inversion symmetry), providing
            symmetry under the three transformations:
                4barz  (4-fold rotation symmetry around z, with inversion)
                3xyz   (3-fold rotation symmetry around xyz-diagonal)
                2z     (2-fold rotation symmetry around z, no inversion)
            """
            self.numSymmetryOperations = 3
            self.symmetryOperations = np.array([
                RotationMatrix(Axis.Z, Angle.PI_HALF, inv=Axis.Z), # "4barz"
                RotationMatrix(Axis.XYZ, Angle.TWO_PI_THIRD),      # "3xyz"
                RotationMatrix(Axis.Z, Angle.PI)])                 # "2z"

        elif (symmetrygroup == PointSymmetryGroup.cubic_m3m):
            """
            Point-symmetry group m3m (cubic, inversion symmetry), providing
            symmetry under the five transformations:
                4z     (4-fold rotation symmetry around z, no inversion)
                3xyz   (3-fold rotation symmetry around xyz-diagonal)
                1bar   (no rotation symmetry, with inversion of all axes)
                2x     (2-fold rotation symmetry around x, no inversion)
                2z     (2-fold rotation symmetry around z, no inversion)
            """
            self.numSymmetryOperations = 5
            self.symmetryOperations = np.array([
                RotationMatrix(Axis.Z, Angle.PI_HALF),              # "4z"
                RotationMatrix(Axis.XYZ, Angle.TWO_PI_THIRD),       # "3xyz"
                RotationMatrix(Axis.XYZ, Angle.ZERO, inv=Axis.XYZ), # "1bar"
                RotationMatrix(Axis.X, Angle.PI),                   # "2x"
                RotationMatrix(Axis.Z, Angle.PI)])                  # "2z"

        elif (symmetrygroup == PointSymmetryGroup.cubic_432):
            """
            Point-symmetry group 432 (cubic, no inversion symmetry), providing
            symmetry under the four transformations:
                4z     (4-fold rotation symmetry around z, no inversion)
                3xyz   (3-fold rotation symmetry around xyz-diagonal)
                2x     (2-fold rotation symmetry around x, no inversion)
                2y     (2-fold rotation symmetry around y, no inversion)
            """
            self.numSymmetryOperations = 4
            self.symmetryOperations = np.array([
                RotationMatrix(Axis.Z, Angle.PI_HALF),              # "4z"
                RotationMatrix(Axis.XYZ, Angle.TWO_PI_THIRD),       # "3xyz"
                RotationMatrix(Axis.X, Angle.PI),                   # "2x"
                RotationMatrix(Axis.Y, Angle.PI)])                  # "2y"

        elif (symmetrygroup == PointSymmetryGroup.cubic_m3):
            """
            Point-symmetry group m3 (cubic, inversion symmetry), providing
            symmetry under the four transformations:
                2z     (2-fold rotation symmetry around z, no inversion)
                3xyz   (3-fold rotation symmetry around xyz-diagonal)
                1bar   (no rotation symmetry, with inversion of all axes)
                2x     (2-fold rotation symmetry around x, no inversion)
            """
            self.numSymmetryOperations = 4
            self.symmetryOperations = np.array([
                RotationMatrix(Axis.Z, Angle.PI),                   # "2z"
                RotationMatrix(Axis.XYZ, Angle.TWO_PI_THIRD),       # "3xyz"
                RotationMatrix(Axis.XYZ, Angle.ZERO, inv=Axis.XYZ), # "1bar"
                RotationMatrix(Axis.X, Angle.PI)])                  # "2x"

        elif (symmetrygroup == PointSymmetryGroup.triclinic_1):
            """
            Point-symmetry group 1 (triclinic, no inversion symmetry),
            not providing symmetry under any transformations except the
            trivial identity:
                1      (no rotation symmetry, no inversion)
            """
            self.numSymmetryOperations = 1
            self.symmetryOperations = np.array([
                RotationMatrix(Axis.NONE, Angle.ZERO)]) # "1"

        elif (symmetrygroup == PointSymmetryGroup.triclinic_1bar):
            """
            Point-symmetry group 1 (triclinic, inversion symmetry), providing
            symmetry under the single transformation:
                1bar   (no rotation symmetry, with inversion of all axes)
            """
            self.numSymmetryOperations = 1
            self.symmetryOperations = np.array([
                RotationMatrix(Axis.XYZ, Angle.ZERO, inv=Axis.XYZ)]) # "1bar"

        elif (symmetrygroup == PointSymmetryGroup.trigonal_3):
            """
            Point-symmetry group 3 (trigonal, no inversion symmetry), providing
            symmetry under the single transformation 3z:
                3z     (3-fold rotation symmetry around z, no inversion)
            """
            self.numSymmetryOperations = 1
            self.symmetryOperations = np.array([
                RotationMatrix(Axis.Z, Angle.TWO_PI_THIRD)])  # "3z" symmetry

        elif (symmetrygroup == PointSymmetryGroup.trigonal_32_2x):
            """
            Point-symmetry group 32(x) (trigonal, no inversion symmetry,
            with 2-fold symmetry for rotation around the x-axis), providing
            symmetry under the two transformations:
                3z     (3-fold rotation symmetry around z, no inversion)
                2x     (2-fold rotation symmetry around x, no inversion)
            """
            self.numSymmetryOperations = 2
            self.symmetryOperations = np.array([
                RotationMatrix(Axis.Z, Angle.TWO_PI_THIRD),   # "3z" symmetry
                RotationMatrix(Axis.X, Angle.PI)])            # "2x" symmetry

        elif (symmetrygroup == PointSymmetryGroup.trigonal_32_2y):
            """
            Point-symmetry group 32(y) (trigonal, no inversion symmetry,
            with 2-fold symmetry for rotation around the y-axis), providing
            symmetry under the two transformations:
                3z     (3-fold rotation symmetry around z, no inversion)
                2y     (2-fold rotation symmetry around y, no inversion)
            """
            self.numSymmetryOperations = 2
            self.symmetryOperations = np.array([
                RotationMatrix(Axis.Z, Angle.TWO_PI_THIRD),   # "3z" symmetry
                RotationMatrix(Axis.Y, Angle.PI)])            # "2y" symmetry

        elif (symmetrygroup == PointSymmetryGroup.hexagonal_6):
            """
            Point-symmetry group 6 (hexagonal, no inversion symmetry),
            providing symmetry under the single transformation 6z:
                6z     (6-fold rotation symmetry around z, no inversion)
            """
            self.numSymmetryOperations = 1
            self.symmetryOperations = np.array([
                RotationMatrix(Axis.Z, Angle.PI_THIRD)])  # "6z" symmetry

        elif (symmetrygroup == PointSymmetryGroup.hexagonal_6bar):
            """
            Point-symmetry group \bar{6} (hexagonal, inversion symmetry),
            providing symmetry under the single transformation 6barz:
                6barz  (6-fold rotation symmetry around z, with inversion)
            """
            self.numSymmetryOperations = 1
            self.symmetryOperations = np.array([
                RotationMatrix(Axis.Z, Angle.PI_THIRD, inv=Axis.Z)]) # "6barz"

        else:
            raise ValueError("Unknown point-symmetry group %s."%symmetrygroup)

        return

    def summarize(self):
        print('='*80)
        print(textwrap.fill("Symmetry operations associated with %s are:"
                            %self.symmetrygroup.value, 80))
        for n in range(len(self.symmetryOperations)):
            print(("\n%d. %s"%(n+1,self.symmetryOperations[n].opdescr)).upper())
            self.symmetryOperations[n].summarize()
        print("\n",end="")
        return

@dataclass
class Tensor:

    version = "1.0"

    def __init__(self, tensortype=TensorType.POLAR, tensorrank=4, init="indices"):
        """
        Instantiate the tensor object.

        Parameters
        ----------
        tensortype : TensorType, optional
            The type of tensor to be instantiated. The two possible choices
            are:
                TensorType.POLAR ("true tensor"), or
                TensorType.AXIAL ("pseudotensor").
                The default is TensorType.POLAR.
        tensorrank : int, optional
            The rank of the tensor. The default is 4.
        init : str, optional
            Initialize the elements of the tensor. The default is "indices".
                "indices" : The elements are initialized as strings containing
                            the indices of the respective elements, without
                            any parentheses or commas. This way, it is a
                            straightforward matter to keep track of the origin
                            of elements once we start analyze, say, the tensor
                            in a rotated or inverted system.
                "none"    : Initialize all elements as empty (length-zero)
                            strings.
        Returns
        -------
        None.
        """
        debug = False
        self.tensortype = tensortype
        self.tensorrank = tensorrank
        self.elements = sp.MutableDenseNDimArray(np.zeros(3*np.ones(tensorrank, dtype=int), dtype=object))
        if debug:
            print("Init: Allocated symbolic rank-%d tensor of size %s."
                  %(len(self.elements.shape), str(self.elements.shape)))
        for index in self.indextuples():
            symstring = self.indexstring(index,lettering=True,strip=True)
            value = sp.symbols(symstring)
            if debug:
                print("Init: Assigning symbol '%s' -> %s"%(symstring, index))
            self.setElement(index, value)
        return

    def indextuples(self):
        indextuples = np.ndarray(shape=(1,self.tensorrank),dtype=int)
        for k in range(self.tensorrank): indextuples[0][k] = 3
        indextuples = np.ndindex(tuple(np.asarray(indextuples[0])))
        return indextuples

    def getElement(self, index):
        value = self.elements[index]
        return value

    def setElement(self, index, value):
        self.elements[index] = value
        return

    def indexstring(self, index, lettering=True, strip=False):
        """
        Generate a string corresponding to the number indices of an element
        of the tensor.
        Examples:
            (0,0,0,0) -> "xxxx"
            (2,0,1,0) -> "zxyx"

        Parameters
        ----------
        index : [int]
            List of integers (with len(index) = {rank of tensor})
            corresponding to an index tuple for the tensor element.
            Notice that the index integers must be listed starting
            from zero, in the regular Python/C/Fortran convention.
        lettering : bool, optional
            If True, use "x","y","z" for indexing; otherwise,
            use "1","2","3". Notice the offset of one in the returned
            string if numbers are used. The default is True.
        strip : bool, optional
            If true, then omit any enclosing "(" or ")". The default is False.

        Returns
        -------
        indexstring : str
            The corresponding index of the element, expressed as a
            regular string.
        """
        offset = ord('x')  # ASCII code for 'x'
        if (not strip):
            indexstring = "("
        else:
            indexstring = ""
        for k in range(self.tensorrank):
            if lettering:  # Use x,y,z for indexing
                indexstring += chr(index[k]+offset)
            else:          # Use 1,2,3 for indexing
                indexstring += str(index[k]+1)
            if (not strip):
                if k < self.tensorrank-1: indexstring += ","
        if (not strip):
            indexstring += ")"
        return indexstring

    def copy(self, source):
        """
        Copy all elements of the source tensor to the present one. This
        procedure assumes that the tensors are of equal rank and of equal
        type (polar or axial). If found not to be so, an error will be issued.

        Parameters
        ----------
        source : Tensor
            The source from which the tensor will be copied.

        Returns
        -------
        None.
        """
        if source.tensortype != self.tensortype:
            raise ValueError('Expecting source of type %s.'%self.tensortype)
        if source.tensorrank != self.tensorrank:
            raise ValueError('Expecting source of rank %d.'%self.tensorrank)
        for index in self.indextuples():
            self.setElement(index, source.getElement(index))
        return

    def rotate(self, rotationmatrix):
        """
        Given a rotation matrix R, either describing a proper or improper
        rotation (with det(R) = 1 or det(R) = -1, respectively), express all
        elements of the current tensor in the rotated coordinate system.

        The identity for the rotation of polar or axial tensors are, for
        example, described in the doctoral thesis of Fredrik Jonsson, "The
        Nonlinear Optics of Magneto-Optic Media" (The Royal Institute of
        Technology, 2000), ISBN 91-7170-575-9; Available via the National
        Library of Sweden (Kungl. Biblioteket) at
        http://urn.kb.se/resolve?urn=urn:nbn:se:kth:diva-2967

        Parameters
        ----------
        rotationmatrix : RotationMatrix
            The rotation matrix, either describing a proper (det(R) = 1) or
            improper (det(R) = -1) rotation.

        Returns
        -------
        rotatedTensor : Tensor
            The corresponding tensor with elements expressed in the rotated
            coordinate system.
        """
        det = rotationmatrix.det()
        rotatedTensor = Tensor(self.tensortype, self.tensorrank, init="none")
        for rotIndex in rotatedTensor.indextuples():
            element = 0
            for index in self.indextuples():
                term = self.getElement(index)
                for n in range(self.tensorrank):
                    term *= rotationmatrix.array[rotIndex[n],index[n]]
                element += term
            if self.tensortype == TensorType.AXIAL:
                if det < 0.0:
                    element = det*element;
            rotatedTensor.setElement(rotIndex, element)
        return rotatedTensor

    def summarize(self, text):
        maxColumnSize = 90
        print('='*80)
        if text != '': print(text)
        print("Tensor is of rank %d and type %s."%(self.tensorrank, str(self.tensortype)))
        print("Tensor consists of %d non-zero elements."%(self.nonzeroElements()))
        print("Tensor elements:")
        column = 0
        for index in self.indextuples():
            block = "%s: %-6s"%(self.indexstring(index), self.getElement(index))
            blocklen = len(block)
            if column + blocklen > maxColumnSize:
                print("\n",end="")
                column = 0
            print(block, end="")
            column += blocklen
        print("\n",end="")
        return

def nonzeroTensorElements(pointsymmetry=PointSymmetryGroup.cubic_23,
            tensortype=TensorType.POLAR,
            tensorrank=4):

    debug = False
    eqs = []
    INDEPENDENT = "independent"    # Used as marker of independent elements
    symmetryoperations = PointSymmetryOperations(pointsymmetry)
    symmetryoperations.summarize()
    tt = Tensor(tensortype, tensorrank, init="indices")

    """
    Build a dictionary containing all elements (also the ones which later on
    might turn out to be zero), for the sake of identification of unaffected
    non-zero elements later on.
    """
    allElements = {}
    for index in tt.indextuples():
        allElements.update({tt.getElement(index): sp.symbols(INDEPENDENT)})

    """
    Build a system of equations by iterating over all symmetry operations
    belonging to the target point-symmetry group, at which the tensor in the
    rotated and/or inverted coordinate system should be equal to the original.
    """
    print(textwrap.fill("Compiling equation system from the %d symmetries of "
          "point-symmetry group %s."%(symmetryoperations.numSymmetryOperations,
                                      pointsymmetry.value),80))
    for n in range(symmetryoperations.numSymmetryOperations):
        so = symmetryoperations.symmetryOperations[n]
        ttr = tt.rotate(so)
        if debug:
            ttr.summarize("Tensor after symmetry operation %s (%s) of %s"
                %(so.oplabel, so.opdescr, symmetryoperations.symmetrygroup))

        """
        Create the system of equations to solve for nonzero and independent
        tensor elements. This forms the core of the automated identification
        of elements in direct inspection. Trivial identities (for example,
        that "zzzz" == "zzzz" under proper rotation around the z-axis, say
        for the trigonal group "3") will hence not be added to the system.
        Therefore, for elements not appearing in the complete system of
        equations (including all restrictions imposed by all symmetry
        operations of the group), these will later on in the extraction of
        nonzero and unique elements hence also be considered as such.
        """
        for index in tt.indextuples():
            eq = tt.getElement(index) - ttr.getElement(index)
            if eq != 0:  # Avoid adding trivial relations to the system
                eqs.append(eq)

    """
    Solve the system of equations for direct inspection and keep only
    non-zero elements for further analysis. At this very stage it would
    be nice to also include removal of linearly dependent (reduntant)
    equations, but it might actually already be done by sympy.solve(eqs)
    before entering the actual solving. (Would be highly interesting to
    know if this is the case.)
    
    Notice that when sorting out non-zero elements, we may still have elements
    which are nonzero but due to symmetry constraints only would yield a
    trivial relation (and hence are omitted from the system of equations).
    Therefore, in order to keep track of these elements later on, we also
    need to sort out the zero-valued elements for later comparison.
    """
    zeroElements = {}
    nonzeroElements = {}
    if len(eqs) == 0:
        print("No linearly independent equations found from system.")
        for key, value in allElements.items():
            nonzeroElements.update({key: value})
    else:
        print("Solving system of %d equations for nonzero elements."%len(eqs))
        sols = sp.solve(eqs)
        if debug: sp.pprint(sols)
        for key, value in sols.items():
            if value != 0:
                nonzeroElements.update({key: value})
            else:
                zeroElements.update({key: value})
        if debug:
            for element, solution in nonzeroElements.items():
                print("%s = %s"%(element, solution))

    """
    Sort out the unique elements from the nonzero set. As the independent
    element, we choose the first one as appearing in the list of index tuples
    for the tensor (typically starting with element 'xxxx'). The unique element
    will be used as a key in the resulting dictionary, with dependent elements
    as values in a list connected to that key.
    """
    uniqueElements = {}
    removed = {}

    for key in nonzeroElements.keys():
        removed[key] = False
    for index, (element, solution) in enumerate(nonzeroElements.items()):
        if (not removed[element]) and (element not in uniqueElements.keys()):
            uniqueElements.update({element: [solution]})
        for testIndex, (testElement, testSolution) in enumerate(nonzeroElements.items()):
            if index < testIndex:
                if (not removed[element]) and (not removed[testElement]):
                    if isinstance(testSolution, sp.Expr):
                        for usol in uniqueElements[element]:
                            for testSymbol in testSolution.free_symbols:
                                if (testSymbol in usol.free_symbols) and (str(testSymbol) != INDEPENDENT):
                                    if debug:
                                        print(textwrap.fill("Found expression '%s' (in "
                                              "expression for unique element '%s = %s') "
                                              "to have a symbol '%s' which also is present "
                                              "in the solution for element '%s = %s'."
                                              %(testSolution, testElement,testSolution,
                                                testSymbol, element, usol), 80));
                                    """
                                    Form a temporary equation and solve with respect to
                                    the dependent expression we found for the previously
                                    enlisted independent element.
                                    """
                                    tsol = sp.solve(testElement - testSolution, testSymbol)
                                    tsol = tsol[0]
                                    if debug:
                                        print("Solving %s = %s with respect to %s, resulting in %s = %s."
                                              %(testElement, testSolution, testSymbol, testSymbol, tsol))
                                    ttsol = usol.subs(testSymbol, tsol)
                                    if debug:
                                        print("Substituting %s into '%s = %s' yields "
                                              "%s = %s."%(tsol, element, usol, element, ttsol))
                                    uniqueElements[element] += [ttsol]
                                    if debug:
                                        print("Added %s as additional dependent equality to "
                                              "list contained by element %s."%(tsol, element))
                                        print("After update, %s = %s."%(element, uniqueElements[element]))
                                        print("Removing the relation %s = %s from "
                                              "the set of independent solutions."
                                              %(testElement, testSolution))
                                    removed[testElement] = True
                                    break
                                break

    """
    To the consolidated set as obtained from the analysis of nonzero and
    independent elements, we should now also add any elements which turned
    out not to be part of the direct inspection, namely any elements which
    were completely unaffected by the set of transformations. Just to mention
    one example, the 3 (trigonal) point-symmetry group only involves rotation
    around the z-axis for the symmetry operation, and hence has its "zzz...z"
    (of arbitrary rank and polar/axial) element completely unaffected.
    """
    for element in allElements.keys():
        if element not in zeroElements.keys():
            if element not in nonzeroElements.keys():
                found = False
                for val in nonzeroElements.values():
                    if element in val.free_symbols:
                        found = True
                if found == False:
                    uniqueElements.update({element: ["independent"]})

    allElements = {}
    for element, value in nonzeroElements.items():
        allElements[element] = True;
        for val in nonzeroElements[element].free_symbols:
            if str(val) != INDEPENDENT:
                allElements[val] = True;
    for element in uniqueElements.keys():
        allElements[element] = True;

    numElements = 3**tensorrank
    numUniqueElements = len(uniqueElements)
    numNonzeroElements = len(allElements)
    print("="*80)
    print(("Point symmetry group %s\n"%pointsymmetry.value).upper())
    print(textwrap.fill("For the %s tensor of rank %d (for which there in a "
        "completely non-symmetrical case would be a maximum of %d elements) "
        "under constraint of point-symmetry group %s, there are %d nonzero "
        "elements, of which %d are independent."
        %(tensortype.value, tensorrank, numElements, pointsymmetry.value,
          numNonzeroElements, numUniqueElements), 80))
    if numUniqueElements > 0:
        print("\nThe %d independent and nonzero tensor elements are as "
              "follows:"%numUniqueElements);
        for index, (key, value) in enumerate(uniqueElements.items()):
            print("    %4d.  %s"%(index+1, key), end="")
            dependentElements = sorted(value, key=lambda x: str(x))
            dependentElements = sorted(dependentElements, key=lambda x: shannon(str(x)))
            for element in dependentElements:
                print(" = %s"%element, end="")
            print("\n", end="")
        print("="*80)
    return

def parseopts(args: List[str]) -> Tuple[str, List[int]]:
    """
    Parse options from command line.

    Parameters
    ----------
    args : List[str]
        The standard argument vector as supplied by the command line options.

    Returns
    -------
    pointsymmetry, tensortype, tensorrank
        pointsymmetry: The point-symmetry group for which we are about to
                       compute nonzero and independent tensor elements.
        tensortype:    Type of tensor (polar or axial)
        tensorrank:    Rank of the tensor
    """
    symmetryOptions = ""
    for s in PointSymmetryGroup:
        symmetryOptions += " "*24
        symmetryOptions += "%-14s- %s"%(s.name, s.value)
        symmetryOptions += "\n"
    typeOptions = ""
    for t in TensorType:
        typeOptions += " "*24
        typeOptions += "%-14s- %s"%(t.name, t.value)
        typeOptions += "\n"
    USAGE = sys.argv[0]+" [OPTIONS]\n"\
            "  -v, --version      List program version and exit.\n"\
            "  -h, --help         List program version and exit.\n"\
            "  -s, --symmetry <s> Specify the point-symmetry group for which\n"\
            "                     the nonzero and independent tensor elements\n"\
            "                     are to be computed. Valid options <s> are:\n"\
            +symmetryOptions+\
            "  -t --type <t>      Specify the tensor type to apply. Valid\n"\
            "                     options <t> are:\n"\
            +typeOptions+\
            "  -r --rank <n>      Specify the rank of the tensor, as an integer <n>."
    pointsymmetry = None
    tensortype = None
    tensorrank = None
    try:
        options, arguments = getopt.getopt(args,
            "vhs:t:r:",
            ["version", "help", "symmetry=", "type=", "rank="])
    except:
        print("Unrecognized option in %s"%args)
        print(USAGE)
        sys.exit()
    for opt, arg in options:
        if opt in ("-v", "--version"):
            print("This is %s v.%s."%(sys.argv[0], Tensor.version))
            sys.exit()
        if opt in ("-h", "--help"):
            print(USAGE)
            sys.exit()
        if opt in ("-s", "--symmetry"): pointsymmetry = PointSymmetryGroup[arg]
        if opt in ("-t", "--type"): tensortype = TensorType[str(arg).upper()]
        if opt in ("-r", "--rank"): tensorrank = int(arg)
    if None in (pointsymmetry, tensortype, tensorrank):
        if (pointsymmetry == None):
            print("Unspecified point symmetry group.")
        if (tensortype == None):
            print("Unspecified tensor type.")
        if (tensorrank == None):
            print("Unspecified tensor rank.")
        print(USAGE)
        sys.exit()
    return pointsymmetry, tensortype, tensorrank

def main() -> None:
    args = sys.argv[1:]
    if len(args) > 0:  # If command-line is supplied
        pointsymmetry, tensortype, tensorrank = parseopts(args)
    else:
        # pointsymmetry=PointSymmetryGroup.cubic_23
        # pointsymmetry = PointSymmetryGroup.trigonal_3
        # pointsymmetry=PointSymmetryGroup.cubic_4bar3m
        # pointsymmetry = PointSymmetryGroup.triclinic_1bar
        # pointsymmetry = PointSymmetryGroup.triclinic_1
        # pointsymmetry = PointSymmetryGroup.tetragonal_422
        # pointsymmetry = PointSymmetryGroup.cubic_m3m
        # pointsymmetry = PointSymmetryGroup.cubic_432
        # pointsymmetry = PointSymmetryGroup.cubic_m3
        # pointsymmetry = PointSymmetryGroup.hexagonal_6 # <= GER FEL ANTAL TERMER!
        # pointsymmetry = PointSymmetryGroup.hexagonal_6bar
        # pointsymmetry = PointSymmetryGroup.trigonal_32_2x
        pointsymmetry = PointSymmetryGroup.trigonal_32_2y
        tensortype = TensorType.POLAR
        tensorrank = 4
    reducedset = nonzeroTensorElements(pointsymmetry, tensortype, tensorrank)
    return

if __name__ == "__main__":
    main()
