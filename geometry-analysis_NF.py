import numpy
import os
import argparse

def calculate_distance(atom1_coord, atom2_coord):
    """
    Calculate the distance between two three-dimensional points.
    """
    x_distance = atom1_coord[0] - atom2_coord[0]
    y_distance = atom1_coord[1] - atom2_coord[1]
    z_distance = atom1_coord[2] - atom2_coord[2]
    bond_length_12 = numpy.sqrt(x_distance ** 2 + y_distance ** 2 + z_distance ** 2)
    return bond_length_12

def bond_check(atom_distance, minimum_length=0, maximum_length=1.5):
    """Check if a distance is a bond based on a minimum and maximum bond length"""

    if atom_distance > minimum_length and atom_distance <= maximum_length:
        return True
    else:
        return False
  

if __name__== "__main__":

    parser = argparse.ArgumentParser(description="This script analyzes a user provided xyz file and outputs all the bonds.")
    parser.add_argument("xyz_file",help="The filepath for the xyz file to analyze")
    args = parser.parse_args()

    file_location = args.xyz_file
    xyz_file = numpy.genfromtxt(fname=file_location, skip_header=2, dtype='unicode')
    symbols = xyz_file[:, 0]
    coordinates = xyz_file[:, 1:]
    coordinates = coordinates.astype(numpy.float)
    num_atoms = len(symbols)
    for num1 in range(0, num_atoms):
        for num2 in range(0, num_atoms):
            if num1 < num2:
                x_distance = coordinates[num1, 0] - coordinates[num2, 0]
                y_distance = coordinates[num1, 1] - coordinates[num2, 1]
                z_distance = coordinates[num1, 2] - coordinates[num2, 2]
                bond_length_12 = numpy.sqrt(x_distance ** 2 + y_distance ** 2 + z_distance ** 2)
                if bond_length_12 > 0 and bond_length_12 <= 1.5:
                    print(F'{symbols[num1]} to {symbols[num2]} : {bond_length_12:.3f}')