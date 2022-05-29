import numpy as np
from asa_dsenum.permutation import gene_trans, get_trans_perms, unique, superperiodic_unique, shin_get_permutation, unique_nov
from reo3_graph import get_weighted_edges_reo3
from rutile_enumerate import get_weighted_edges_rutile
from rutile_oct_enumrate import get_weighted_edges_rutile_with_oct


def get_weighted_edges(st, mode='rutile', eps=0.001, valence=2):
    if mode == 'rutile':
        return get_weighted_edges_rutile(st, eps=eps, valence=valence)
    elif mode == 'rutile_with_oct':
        return get_weighted_edges_rutile_with_oct(st, eps=eps, valence=valence)
    elif mode == 'reo3':
        return get_weighted_edges_reo3(st, eps=eps, valence=valence)
    else:
        raise KeyError


def get_all_vectors(gs_all, supercell):
    vectors = []
    for graph in gs_all:
        vectors.append(graph_to_vector(graph, supercell))

    return vectors


def graph_to_vector(graph, supercell):

    vertexs = get_vertexs(graph)
    onehot_vector = ''
    for i, site in enumerate(supercell.sites):
        if site.specie.name == 'O':
            if i in vertexs:
                onehot_vector += '1'
            else:
                onehot_vector += '0'

    return onehot_vector


def get_vertexs(graph):

    vertex_group = set()
    for i, j in graph:
        vertex_group.add(i)
        vertex_group.add(j)

    return vertex_group


def o_sublattice_from_rutile(structure, mode='rutile'):
    atomic_nums = list(structure.atomic_numbers)
    if mode == 'rutile':
        cation_nums = int(len(atomic_nums)/3)
    elif mode == 'reo3':
        cation_nums = int(len(atomic_nums)/4)
    else:
        raise KeyError
    lattices = structure.lattice.matrix
    fractional_coords = structure.frac_coords[cation_nums:]

    return [lattices, fractional_coords, atomic_nums[cation_nums:]]


def get_permutations(supercell, parent_sym, hnf, mode='rutile'):

    o_sublattice = o_sublattice_from_rutile(supercell, mode=mode)

    dic_zahyou = dict()  # 分率座標の行列表示
    for i in range(len(o_sublattice[1])):
        dic_zahyou[i] = np.asarray(o_sublattice[1][i])

    permutations = shin_get_permutation(
        o_sublattice, parent_sym, hnf)
    set_of_translations = gene_trans(hnf)
    SuperPeriodicPermutations = get_trans_perms(
        dic_zahyou, set_of_translations)

    return permutations, SuperPeriodicPermutations


def graph_unique(graph_all, supercell, parent_sym, hnf, mode='rutile'):

    permutations, spp = get_permutations(
        supercell, parent_sym, hnf, mode=mode)

    struct_candidates = get_all_vectors(graph_all, supercell)

    # unique
    structures = unique(permutations, struct_candidates)
    # structures = unique_nov(struct_candidates, permutations)

    # superperiodic deletion
    non_sp_structures = superperiodic_unique(spp, structures)

    return structures, non_sp_structures


def graph_unique_nov(graph_all, supercell, parent_sym, hnf, mode='rutile'):

    permutations, spp = get_permutations(
        supercell, parent_sym, hnf, mode=mode)

    struct_candidates = get_all_vectors(graph_all, supercell)
    print(f'num of candidate = {len(struct_candidates)}')

    # unique
    # structures = unique(permutations, struct_candidates)
    structures = unique_nov(permutations, struct_candidates)

    # superperiodic deletion
    non_sp_structures = superperiodic_unique(spp, structures)

    return structures, non_sp_structures
