"""
molssi-geometry-analysis.py
This module contains the outcome of the MolSSI geometry analysis project
Requires one command-line argument specifying an xyz file to analyze
"""

import numpy, os, sys

def open_xyz(filename):
    '''
    Open the specified xyz file
    Parameters: The filename of the xyz file
    Returns: The symbols (list); the atomic coordinates (numpy array)
    '''
    # Get raw data from xyz file
    xyz_file = numpy.genfromtxt(fname=filename,skip_header=2,dtype='unicode')

    # Extract symbols and coordinates from raw data
    symbols = xyz_file[:,0]
    coordinates = xyz_file[:,1:]
    coordinates = coordinates.astype(numpy.float)
    return symbols, coordinates

def bond_check(atom1_coord,atom2_coord,min_length=0,max_length=1.5):
    '''
    Checks whether the distance between two atoms is greater than min_length
    and less than max_length
    Required Parameters: atom1_coord, atom2_coord
    Optional Parameters: min_length (default 0 Angstrom),
                     max_length (default 1.5 Angstrom)
    Returns: True or False
    '''
    distance = calculate_distance(atom1_coord,atom2_coord)
    # Confirm that the distance is non-negative
    if distance < 0:
        raise ValueError(f"negative distance detected ({distance})")
    if distance > min_length and distance < max_length:
        return True
    else:
        return False

def calculate_distance(atom1_coord,atom2_coord):
    """
    Computes the distance between the two atoms with coordinates specified by atom1_coord and atom2_coord
    Parameters: atom1_coord, atom2_coord
    Returns: the distance between the atoms
    """
    # Check that coordinates have the correct max_length
    if (len(atom1_coord) != 3 or len(atom2_coord) != 3):
        raise ValueError("Atomic coordinates with incorrect length detected")
    x_distance = atom2_coord[0]-atom1_coord[0]
    y_distance = atom2_coord[1]-atom1_coord[1]
    z_distance = atom2_coord[2]-atom1_coord[2]
    r_distance = numpy.sqrt(x_distance**2 + y_distance**2 + z_distance**2)
    return r_distance

#def geometry_analysis():
if __name__ == "__main__":
    # Check that a filename has been provided as an argument
    if len(sys.argv) < 2:
        raise NameError("Missing xyz file argument. Please specify an xyz file to be analyzed!")
    file_location = sys.argv[1]
    #Alternative approach: use try/except
    #try:
    #    file_location = sys.argv[1]
    #except:
    #    raise NameError("Unable to interpret xyz file location. Please specify an xyz file to be analyzed!")

    # Alternative xyz file acquisition (hard-coded)
    #file_location = os.path.join('Desktop','cms-workshop','data','data','water.xyz')

    symbols, coordinates = open_xyz(file_location)
    num_atoms = len(symbols)

    # Compute the matrix of distances
    for i in range(num_atoms):
        for j in range(i,num_atoms):
            distance = calculate_distance(coordinates[i],coordinates[j])
            if bond_check(coordinates[i],coordinates[j]):
                print(F'{symbols[i]} to {symbols[j]} : {distance:.3f}')

#if __name__ == "__main__":
#    geometry_analysis()
