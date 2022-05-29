import generate_HNF
import generate_superlattice
import permutation
import numpy as np
import spglib
from copy import deepcopy
import time

rutile = ([(4, 0, 0), (0, 4, 0), (0, 0, 3)],
          [(0, 0, 0), (0.5, 0.5, 0.5), (0.3, 0.3, 0.0),
           (0.7, 0.7, 0.0), (0.2, 0.8, 0.5), (0.8, 0.2, 0.5)],
          [50, 50, 8, 8, 8, 8])


def rutile_enumration(base_structure, index, sublattice=False):

    # 処理前の時刻
    t1 = time.time()
    parent_sym = spglib.get_symmetry(base_structure)
    HNF_list = generate_HNF.generate_all_superlattices(index)
    reduced_HNF = generate_HNF.reduce_HNF_list_by_parent_lattice_symmetry(
        HNF_list, parent_sym["rotations"])
    for number, matrix in enumerate(reduced_HNF):
        ds_list = list()

        parent_lattice = generate_superlattice.get_superlattice(
            rutile, matrix, index)
        o_sublattice = generate_superlattice.get_o_superlattice(
            rutile, matrix, index)
        set_of_translations = permutation.gene_trans(matrix)
        dic_zahyou = dict()
        for i in range(len(o_sublattice[1])):
            dic_zahyou[i] = o_sublattice[1][i]
        set_of_transtikans = permutation.get_trans_perms(
            dic_zahyou, set_of_translations)

        goalen = permutation.shin_get_permutation(
            o_sublattice, parent_sym, matrix)
        for fuga in range(1, 4*index):
            pattern = permutation.make_candidate(o_sublattice, fuga)
            uniqued_pattern = permutation.unique(goalen, pattern)
            unique_structures = permutation.superperiodic_unique(
                set_of_transtikans, uniqued_pattern)
            ds_list.extend(unique_structures)
        number += len(ds_list)

    print('見つかったのは{}個の構造パターン'.format(number))

    # 処理後の時刻
    t2 = time.time()

    # 経過時間を表示
    elapsed_time = t2-t1
    print("かかった時間は", elapsed_time, "秒")


if __name__ == '__main__':
    rutile_enumration(rutile, 2)
