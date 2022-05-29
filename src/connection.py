import os
from copy import deepcopy
import time
from pymatgen.io.cif import CifWriter
from pymatgen.io.cif import CifParser
from pymatgen.core import Specie, DummySpecie
from pymatgen.core import Lattice, Structure
from asa_dsenum import generate_HNF, generate_superlattice, permutation
import numpy as np
import spglib
import pymatgen

mapping_color_species = [DummySpecie('X'), Specie('O'), Specie('Sn')]

fir_vecs = [(0.4, 0.4, 0),
            (-0.1, 0.5, 0.5), (0.5, -0.1, 0.5), (-0.1, -0.5, 0.5), (-0.5, -0.1, 0.5),
            (-0.1, 0.5, -0.5), (0.5, -0.1, -
                                0.5), (-0.1, -0.5, -0.5), (-0.5, -0.1, -0.5),
            (0.1, -0.5, 0.5), (-0.5, -0.1, 0.5), (0.1, 0.5, 0.5), (0.5, 0.1, 0.5),
            (0.1, -0.5, -0.5), (-0.5, -0.1, -
                                0.5), (0.1, 0.5, -0.5), (0.5, 0.1, -0.5),
            (0.4, -0.4, 0)
            ]


def generate_void_list(mate):
    """generate void list from mate"""
    void_list = list()
    for j in range(len(mate)):
        if int(mate[j]) == 0:
            void_list.append(j)
    return void_list


def coordinate(i, parent_lattice):
    """transform fracitional coordinate to cartersian coordinate"""
    A1_vecs = np.array(parent_lattice[0][0])
    A2_vecs = np.array(parent_lattice[0][1])
    A3_vecs = np.array(parent_lattice[0][2])
    new_coordinate = A1_vecs * parent_lattice[1][i][0] + A2_vecs * \
        parent_lattice[1][i][1] + A3_vecs * parent_lattice[1][i][2]
    return new_coordinate


def calc_distance(coordinate1, coordinate2):
    distance = (coordinate1[0]-coordinate2[0])**2 + (coordinate1[1] -
                                                     coordinate2[1])**2 + (coordinate1[2] - coordinate2[2])**2
    return distance


def distance_check(coordinate1, coordinate2, shortest_length=2.7):
    """check whether if distance between coorinate1 and coordinate2 above shortest distance
    coordinate is a numpy array"""
    counter = False
    if calc_distance(coordinate1, coordinate2) < shortest_length:
        counter = True
    return counter


def connected_check(o_sublattice, mate):
    void_list = generate_void_list(mate)
    unconnected_void_list = list()
    for void_num in void_list:
        break_swich = False
        void1_coordinate = coordinate(void_num, o_sublattice)
        for gensi_num in void_list:
            if void_num == gensi_num:
                break
            void2_coordinate = coordinate(gensi_num, o_sublattice)
            if connected(void1_coordinate, void2_coordinate):
                break_swich = True
                break
        if break_swich:
            print("free")
        else:
            unconnected_void_list.append(void_num)
        return unconnected_void_list
