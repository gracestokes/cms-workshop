"""
This module has functions associated with analyzing the geometry of a molecule. It does not work properly yet. I use geometry-analysis_NF.py instead.
"""

import numpy
import os
import argparse

def calculate_distance(coordinates1,coordinates2):
    x_distance=coordinates1[0]-coordinates2[0]
    y_distance=coordinates1[1]-coordinates2[1]
    z_distance=coordinates1[2]-coordinates2[2]
    distance_12 = numpy.sqrt(x_distance**2+y_distance**2+z_distance**2)
    return distance_12

def bond_check(atom_distance, minimum_length=0, maximum_length=1.5):
    """
    Check if a distance is a bond based on a minimum and maximum bond length
    """
    if atom_distance > minimum_distance and atom_distance <= maximum_length:
        return True
    else:
        return False
    #True and False have a defined meaning in Python (boolean operators)

if __name__== "__main__":
    
    parser = argparse.ArgumentParser(description="This scrip analyzes a user given xyz file and outputs the length of the bonds."
    parser.add_argument("xyz_file", help="The filepath for the xyz file to analyze")
    args=parser.parse_args()                               

file_location=args.xyz_file
xyz_file=numpy.genfromtxt(fname=file_location, skip_header=2, dtype='unicode')
symbols = xyz_file[:,0]
coordinates = xyz_file[:,1:]
coordinates = coordinates.astype(numpy.float)
num_atoms = len(symbols)
for atom1 in range(0,num_atoms):
    for atom2 in range(0,num_atoms):
        if atom1 < atom2:
            x_distance=coordinates1[0]-coordinates2[0]
            y_distance=coordinates1[1]-coordinates2[1]
            z_distance=coordinates1[2]-coordinates2[2]
            distance_12 = numpy.sqrt(x_distance**2+y_distance**2+z_distance**2)
            #this calls the function below and you need to have already SHIFT-ENTERED the cell below for this cell to work
            if distance_12 is True:
                #is True is optional because default is True
                print(F'{symbols[atom1]} to {symbols[atom2]}:{distance_12:.3f}')        