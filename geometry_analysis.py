"""
geometry_analysis.py
This module contains the geometry analysis project.

Author: Paul A. Craig
"""

# Format: docstring, imports, functions, code
import numpy
import os
import sys

def open_xyz(file_location):
    """
    Opens and processes a standard format xyz file and returns two outputs: the symbols and a numpy array of coordinates
    """
    xyz_file = numpy.genfromtxt(fname = file_location, skip_header = 2, dtype = 'unicode')
    # you must use skip header to get rid of lines that don't contain the same number of columns

    # Separate data types
    symbols = xyz_file[:,0] #symbols are already strings
    coordinates = xyz_file[:,1:]
    coordinates = coordinates.astype(numpy.float)
    return symbols, coordinates

def bond_check(distance, minimum_length=0, maximum_length=1.5):
    """
    This function determines if a distance is between the specified minimum and maximum values.
    The specified minimum is 0 Angstroms and the default maximum value is 1.5 Angstroms.It returns true or false
    Input: distance, minimum_length, maximum_length
    Return: Boolean True or False
    """
    if distance > minimum_length and distance <= maximum_length:
        return True
    else:
        return False
    if distance<0:
        raise ValueError(F"Bond Check has detected a negative distance: {distance} Angstroms. Check your input.")

def calculate_distance(atom1_coord, atom2_coord):
    """
    This function takes to coordinates of two atoms and calculates the distance between them.
    Inputs: atom1_coord, atom2_coord
    Return: distance
    """
    x_distance = atom1_coord[0] - atom2_coord[0]
    y_distance = atom1_coord[1] - atom2_coord[1]
    z_distance = atom1_coord[2] - atom2_coord[2]
    distance = numpy.sqrt(x_distance**2 + y_distance**2 + z_distance**2)
    return distance # must be something calculated in the body of code

#file_location = os.path.join('data', 'water.xyz')

if __name__ == "__main__":
    print(F"Running {sys.argv[0]}.")
    if len(sys.argv) < 2:
        raise NameError("Incorrect input! Please specify an xyz file to be analyzed.")
    file_location = sys.argv[1] # You need to provide a file name to be analyzed.

    symbols, coord = open_xyz(file_location)
    num_atoms = len(symbols)
    for num1 in range(0,num_atoms):
        for num2 in range(0,num_atoms):
            if num1<num2:
                bond_length_12 = calculate_distance(coord[num1], coord[num2])
                if bond_check(bond_length_12) is True:
                    print(F'{symbols[num1]} to {symbols[num2]} : {bond_length_12:.3f}')
