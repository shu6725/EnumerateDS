# generating superlattice
import numpy as np
import generate_HNF
import spglib
from copy import deepcopy


def get_superlattice_vectors(o_sublattice, HNF):
    """

    generate superlattice-vector from o_sublattice and HNF


    Parameters

    ----------

    o_sublattice : list[]

    HNF : list[list[]]3 * 3 integer matrix

    Returns

    -------

    list_HNF: list of matrices, each element is Hermite normal form

    """
    lattice = np.asarray(o_sublattice[0])
    lattice = lattice.dot(HNF)
    lattice = lattice.T
    lattice = lattice.tolist()
    return lattice


def shift(coordinate, a, i):
    original = deepcopy(coordinate)
    for j in range(1, a[i]):
        vec = np.zeros(3)
        vec[i] = j
        t = deepcopy(original)
        for k in range(len(original)):
            t[k] += vec
        coordinate = np.concatenate([coordinate, t])
    return coordinate


def create_superlattice(superlattice, HNF):
    a = np.sum(HNF, axis=1)
    for i in range(3):
        superlattice = shift(superlattice, a, i)
    return superlattice

# 　superlattice を生成


def making_superlattice(HNF):
    Sn_sublattice = [(0, 0, 0), (0.5, 0.5, 0.5)]
    O_sublattice = [(0.3, 0.3, 0.0), (0.7, 0.7, 0.0),
                    (0.2, 0.8, 0.5), (0.8, 0.2, 0.5)]
    Sn_superlattice = create_superlattice(Sn_sublattice, HNF)
    O_superlattice = create_superlattice(O_sublattice, HNF)
    superlattice = np.concatenate([Sn_superlattice, O_superlattice])

    copy_superlattice = deepcopy(superlattice)
    renew = np.linalg.inv(HNF)  # 逆行列 新しいsuperlatticevector で分立座標にあらわせるように　
    renew = renew.T
    for i in range(len(superlattice)):
        copy_superlattice[i] = copy_superlattice[i].dot(renew)

    # 消すべき行を記録
    delete_list = []
    for garbage in range(len(copy_superlattice)):
        if 0 <= copy_superlattice[garbage][0] < 0.99 and 0 <= copy_superlattice[garbage][1] < 0.99 and 0 <= copy_superlattice[garbage][2] < 0.99:
            continue
        else:
            delete_list.append(garbage)

    new_zahyou = np.delete(copy_superlattice, delete_list, 0)
    return new_zahyou


def making_o_subsuperlattice(HNF):
    O_sublattice = [(0.3, 0.3, 0.0), (0.7, 0.7, 0.0),
                    (0.2, 0.8, 0.5), (0.8, 0.2, 0.5)]
    O_superlattice = create_superlattice(O_sublattice, HNF)
    superlattice = O_superlattice

    copy_superlattice = deepcopy(superlattice)
    renew = np.linalg.inv(HNF)  # 逆行列 新しいsuperlatticevector で分立座標にあらわせるように　
    renew = renew.T
    for i in range(len(superlattice)):
        copy_superlattice[i] = copy_superlattice[i].dot(renew)

    # 消すべき行を記録
    delete_list = []
    for garbage in range(len(copy_superlattice)):
        if 0 <= copy_superlattice[garbage][0] < 0.99 and 0 <= copy_superlattice[garbage][1] < 0.99 and 0 <= copy_superlattice[garbage][2] < 0.99:
            continue
        else:
            delete_list.append(garbage)

    new_zahyou = np.delete(copy_superlattice, delete_list, 0)
    return new_zahyou


def get_parent_gensi(index):
    empty = []
    for i in range(index*2):
        empty.append(50)
    for i in range(index*4):
        empty.append(8)
    return empty


def get_o_sub_gensi(index):
    empty = []
    for j in range(index*4):
        empty.append(8)
    return empty


def get_superlattice(parent_lattice, HNF, index):
    """
        spglibに適合するフォーマットのsuperlatticeを得る
    """
    superlattice = [[], [], []]  # a は　新しいsuperlattice
    superlattice[0] = get_superlattice_vectors(parent_lattice, HNF)
    superlattice[1] = making_superlattice(HNF)
    superlattice[1] = np.round(superlattice[1], 4)
    superlattice[1] = superlattice[1].tolist()
    superlattice[2] = get_parent_gensi(index)
    superlattice = tuple(superlattice)
    return superlattice


def get_o_superlattice(parent_lattice, HNF, index):
    """
        spglibに適合するフォーマットの酸素superlatticeを得る
    """
    superlattice = [[], [], []]  # a は　新しいsuperlattice
    superlattice[0] = get_superlattice_vectors(parent_lattice, HNF)
    superlattice[1] = making_o_subsuperlattice(HNF)
    superlattice[1] = np.round(superlattice[1], 4)
    superlattice[1] = superlattice[1].tolist()
    superlattice[2] = get_o_sub_gensi(index)
    superlattice = tuple(superlattice)
    return superlattice
