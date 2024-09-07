import os
import shutil
from ase.io import read, write
import numpy as np

def count_oxygen_neighbors(poscar_path, metal_atoms,layer_thickness):
    # Read the POSCAR file
    structure = read(poscar_path)
   # print(len(structure))
    # Determine the z-coordinate range for the topmost layer
    z_coordinates = [atom.position[2] for atom in structure if atom.symbol in metal_atoms]
    max_z = max(z_coordinates)
    min_z_for_top_layer = max_z - layer_thickness
    print(max_z)
    print(min_z_for_top_layer)
    # Identify the topmost metal atoms
   # topmost_metal_atoms = [atom.index for atom in structure if atom.symbol in metal_atoms and atom.position[2] >= min_z_for_top_layer]
    while True:
        topmost_metal_atoms = [atom.index for atom in structure if atom.symbol in metal_atoms and atom.position[2] >= min_z_for_top_layer]
        print(len(topmost_metal_atoms))
        print(topmost_metal_atoms)
        if len(topmost_metal_atoms) < 6:
            print("the number of metal atoms is lower than 6")
        # Increase min_z_for_top_layer and continue the loop
            min_z_for_top_layer -=0.1
        elif len(topmost_metal_atoms) > 6:
            print("the number of metal atoms is higher than 6")
            min_z_for_top_layer = min_z_for_top_layer+0.05
        # Break the loop if the target number of atoms is reached
        elif len(topmost_metal_atoms) == 6:
            break
    top_metal_atom_z =[structure[index].position[2] for index in topmost_metal_atoms]

    indices_of_min_values = sorted(range(len(top_metal_atom_z)), key=lambda i: top_metal_atom_z[i])[:3]

    # Get the corresponding indices in topmost_metal_atoms
    top_three_indices = [topmost_metal_atoms[i] for i in indices_of_min_values]

    # Initialize count dictionary for topmost metal atoms

    unique_metal_types = set()
    unique_min_indices = []
    for index in top_three_indices:
        metal_type = structure[index].symbol
        if metal_type not in unique_metal_types:
            unique_metal_types.add(metal_type)
            unique_min_indices.append(index)

    return unique_min_indices


def add_O(poscar_path, index, directory):
    # Read the POSCAR file
    slab = read(poscar_path)
    total_atom_number = slab.get_global_number_of_atoms()
    metal_ads_site_pos = slab.positions[index]
    slab.append('O')  # add the first oxygen atom
    slab.positions[total_atom_number] = metal_ads_site_pos + [0, 0, 1.8]
    # Save in the specified directory
    slab.write(f'{directory}/POSCAR_O2.vasp')

def add_OH(poscar_path, index, directory):
    # Read the POSCAR file
    slab = read(poscar_path)
    total_atom_number = slab.get_global_number_of_atoms()
    metal_ads_site_pos = slab.positions[index]
    slab.append('O')  # add the first oxygen atom
    slab.positions[total_atom_number] = metal_ads_site_pos + [0, 0, 1.8]
    slab.append('H')  # add the second oxygen atom
    slab.positions[total_atom_number + 1] = slab.positions[total_atom_number] + [0.06, 0.9, 0.5]
    #second_oxygen_pos = metal_ads_site_pos +[0, 0, 1.8]+ [0, 0, 1.246]
    #hydrogen_pos = second_oxygen_pos + [0, 0, 1.05]
    #angle = 109.5
    #rotation_matrix = np.array([[np.cos(np.radians(angle)), -np.sin(np.radians(angle)), 0],
     #                           [np.sin(np.radians(angle)), np.cos(np.radians(angle)), 0],
      #                          [0, 0, 1]])

    # Apply the rotation matrix to the hydrogen position
    #rotated_hydrogen_pos = inp.dot(rotation_matrix, hydrogen_pos - second_oxygen_pos) + second_oxygen_pos
    #slab.append('H')
    # Add the rotated hydrogen atom to the slab
    #slab.append(Atom('H', position=rotated_hydrogen_pos))
    #slab.positions[total_atom_number + 2] = rotated_hydrogen_pos

    # Save in the specified directory
    slab.write(f'{directory}/POSCAR_OH.vasp')


def add_O2(poscar_path, index, directory):
    # Read the POSCAR file
    slab = read(poscar_path)
    total_atom_number = slab.get_global_number_of_atoms()
    metal_ads_site_pos = slab.positions[index]
    slab.append('O')  # add the first oxygen atom
    slab.positions[total_atom_number] = metal_ads_site_pos + [0, 0, 1.8]
    slab.append('O')  # add the second oxygen atom
    slab.positions[total_atom_number + 1] = slab.positions[total_atom_number] + [0, 0, 1.246]

    # Save in the specified directory
    slab.write(f'{directory}/POSCAR_O2.vasp')

def add_OOH(poscar_path, index, directory):
    # Read the POSCAR file
    slab = read(poscar_path)
    total_atom_number = slab.get_global_number_of_atoms()
    metal_ads_site_pos = slab.positions[index]
    slab.append('O')  # add the first oxygen atom
    slab.positions[total_atom_number] = metal_ads_site_pos + [0, 0, 1.8]
    slab.append('O')  # add the second oxygen atom
    slab.positions[total_atom_number + 1] = slab.positions[total_atom_number] + [0.8, 0, 0.7]
    #second_oxygen_pos = metal_ads_site_pos +[0, 0, 1.8]+ [0, 0, 1.246]
    slab.append('H')  # add the hydrogen atom
    slab.positions[total_atom_number + 2] = slab.positions[total_atom_number + 1] + [-0.3, 0, 0.7]
    slab.write(f'{directory}/POSCAR')

   



def create_directory_and_save_files(poscar_path, unique_min_indices):
    #upper_level_path = os.path.join(os.path.dirname(poscar_path), "..")  # One level up from the POSCAR file's directory
    count = 0
    for index in unique_min_indices:
        for add_atom in ["O","OH","OO","HOO"]:
            count += 1
            metal_type = read(poscar_path)[index].symbol
            directory = os.path.join(f"{count}_{metal_type}_{index}_{add_atom}")
                # Create directory if it doesn't exist
            if not os.path.exists(directory):
                os.makedirs(directory)
                if add_atom == "O":
                    # Add O2 and OOH and save in the respective directories
                    add_O(poscar_path, index, directory)
                elif add_atom == "OH":
                    add_OH(poscar_path, index, directory)

                elif add_atom == "OO":
                    add_O2(poscar_path, index, directory)

                elif add_atom == "HOO":
                    add_OOH(poscar_path, index, directory)

metal_atoms = ['Fe', 'Co', 'Ni', 'Zn', 'Sn', 'Sb', 'Bi']  # Replace with your list of metal atom symbols
threshold_distance = 2.5  # Replace with your threshold distance in Angstrom
layer_thickness = 0.5  # Thickness of the topmost layer in Angstrom
unique_min_indices = count_oxygen_neighbors('./CONTCAR', metal_atoms,layer_thickness)
create_directory_and_save_files('./CONTCAR', unique_min_indices)

