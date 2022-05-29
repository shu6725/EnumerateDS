# generate permutation
from copy import deepcopy
import numpy as np
import spglib
from generate_superlattice import get_o_superlattice, get_superlattice


def is_unimodular(M: np.ndarray) -> bool:

    if np.abs(np.around(np.linalg.det(M))) == 1:

        return True

    else:

        return False


def cast_integer_matrix(arr: np.ndarray) -> np.ndarray:

    arr_int = np.around(arr).astype(np.int)

    return arr_int


def is_integer_matrix(matrix):
    for i in range(0, 3):
        for j in range(0, 3):
            if not matrix[i][j].is_integer():
                return False
    return True


def is_site_same(pos1, pos2, ips=3):

    determinant = ((pos1[0]-pos2[0])**2 + (pos1[1]-pos2[1])
                   ** 2 + (pos1[2]-pos2[2])**2) ** (1/2)

    if determinant < 10 ** (-ips):
        return True
    else:
        return False


def get_superlattice_symmetry(parent_sym, HNF):
    H_inv = np.linalg.inv(HNF)
    Subspacegroup = dict()
    Subspacegroup["rotations"] = list()
    Subspacegroup["translations"] = list()
    for i in range(len(parent_sym["rotations"])):
        check = H_inv.dot(parent_sym["rotations"][i].dot(HNF))
        if is_integer_matrix(check):
            t = parent_sym["translations"][i]
            t_neo = np.dot(H_inv, t)
            Subspacegroup["rotations"].append(check)
            Subspacegroup["translations"].append(t_neo)
    return Subspacegroup


def generate_permutation_dict(rotation, translation, dic_zahyou):
    """
    対象操作について置換の移り先を示す辞書を作成

    parameters

    dic_zahyou:元の座標セット

    retruns

    """

    dic2 = deepcopy(dic_zahyou)
    for l in range(len(dic_zahyou)):
        dic2[l] = (rotation.dot(dic_zahyou[l]) + translation)
        dic2[l] = np.round(dic2[l], 4) % 1

    permutation_dict = dict()
    for l in range(len(dic_zahyou)):
        for j in range(len(dic_zahyou)):
            if is_site_same(dic2[l], dic_zahyou[j]):
                permutation_dict[l] = j
    return permutation_dict


def shin_get_permutation(o_super_sublattice, parent_sym, HNF):
    dic_zahyou = dict()  # 分率座標の行列表示
    for i in range(len(o_super_sublattice[1])):
        dic_zahyou[i] = np.asarray(o_super_sublattice[1][i])
    trans = gene_trans(HNF)
    trans_tikans = get_trans_perms(dic_zahyou, trans)
    super_point_group = get_superlattice_symmetry(parent_sym, HNF)

    def generate_super_perm(dic_zahyou, super_point_group):
        # 各座標をnumpy変換してgenerate_permutation_dictを使う
        jun = dict()
        for i in range(len(super_point_group["translations"])):
            jun[i] = generate_permutation_dict(super_point_group["rotations"][i],
                                               super_point_group["translations"][i], dic_zahyou)
        return jun

    rot_permutations = generate_super_perm(dic_zahyou, super_point_group)

    goal_perms = generate_abs_permuatation(rot_permutations, trans_tikans)

    return goal_perms


def generate_super_perm(dic_zahyou, super_point_group):
    # 各座標をnumpyに変換して、generate_permutation_dictを使う
    jun = dict()
    for i in range(len(super_point_group["translations"])):
        jun[i] = generate_permutation_dict(super_point_group["rotations"][i],
                                           super_point_group["translations"][i], dic_zahyou)
    return jun


def gene_trans(HNF):
    """
    HNFから並進操作の一覧を作成

    parameters:HNF

    retruns

    """
    vectors = np.zeros(3)
    H_inv = np.linalg.inv(HNF)
    hantei = True
    for dimention in range(3):
        cp_vecs = deepcopy(vectors)
        if HNF[dimention][dimention] > 1:
            for king in range(1, HNF[dimention][dimention]):
                cp_vecs2 = deepcopy(cp_vecs)
                if hantei:
                    cp_vecs2[dimention] += king
                else:
                    for num_vecs in range(cp_vecs2.shape[0]):
                        cp_vecs2[num_vecs][dimention] += king
                vectors = np.vstack((vectors, cp_vecs2))
            hantei = False
    for hoge in range(len(vectors)):
        vectors[hoge] = np.dot(H_inv, vectors[hoge])
    return vectors


def get_trans_permuation(dic_zahyou, translation):
    """
    並進操作について置換の移り先を示す辞書を作成

    parameters

    retruns

    """
    tin = dict()
    dic2 = deepcopy(dic_zahyou)
    for l in range(len(dic_zahyou)):
        dic2[l] = (dic_zahyou[l] + translation)
        dic2[l] = np.round(dic2[l], 3) % 1
    for l in range(len(dic_zahyou)):
        for j in range(len(dic_zahyou)):
            if is_site_same(dic2[l], dic_zahyou[j]):
                tin[l] = j
    return tin


def get_trans_perms(dic_zahyou, translations):
    """
    並進操作について置換の移り先を示す辞書の辞書を作成

    parameters

    retruns

    """
    trans_perms = dict()
    for i in range(len(translations)):
        trans_perms[i] = get_trans_permuation(dic_zahyou, translations[i])
    return trans_perms


def get_combinations(rot_tikan, trans_tikan):
    """
    回転操作と並進操作の置換を組合す

    parameters

    retruns

    """
    ans = dict()
    for i in range(len(rot_tikan)):
        ans[i] = trans_tikan[rot_tikan[i]]
    return ans


def superperiodic_unique(trans_tikans, omomi4):
    pattern_copy = deepcopy(omomi4)
    lis = set()

    for candidate in range(len(omomi4)):
        if omomi4[candidate] in pattern_copy:
            # ある元の構造が生き残ってるかチェック
            id_superperiodic = False

            for permutation in range(1, len(trans_tikans)):
                # 恒等操作以外の並進操作について回す
                d = dict()
                sta = ""
                for k in range(len(trans_tikans[0])):
                    d[trans_tikans[permutation][k]] = omomi4[candidate][k]
                for r in range(len(trans_tikans[0])):
                    sta += d[r]    # 置換によって生成した配列
                # 作った配列と元の配列が一致、すなわちsuper_periodicやったら
                if int(sta) == int(omomi4[candidate]):
                    id_superperiodic = True
                else:  # 操作によって元とは別の配列になったかどうか
                    if sta in pattern_copy:  # まだ省いていない構造やったら
                        pattern_copy.remove(sta)

            if not id_superperiodic:
                lis.add(omomi4[candidate])

    lis = list(lis)
    return lis


def generate_abs_permuatation(parent_lattice_jun, trans_perms):
    """
    回転操作と並進操作の置換を組合せた究極の辞書を作成

    parameters

    retruns

    """
    tikan_list = list()
    for i in parent_lattice_jun:
        for j in trans_perms:
            tikan = get_combinations(parent_lattice_jun[i], trans_perms[j])
            tikan_list.append(tikan)
    return tikan_list


# 各座標をnumpyに変換して、generate_permutation_dictを使う
def generate_per(superlattice, o_sublattice):
    parent_sym = spglib.get_symmetry(superlattice)
    tin = dict()
    arr = np.asarray(o_sublattice[1])  # lattcie no 行列表示
    dic = dict()  # 分率座標の行列表示
    for i in range(len(o_sublattice[1])):
        dic[i] = np.asarray(o_sublattice[1][i])

    jun = dict()
    for i in range(len(parent_sym["translations"])):
        jun[i] = generate_permutation_dict(parent_sym["rotations"][i],
                                           parent_sym["translations"][i], dic)
    return jun

# 置換の積の式を作る [{0, 1, 2, 3}, {4, 5, 6, 7}]


def make_candidate(o_sublattice, atom_num):
    pattern = []
    for i in range(2**len(o_sublattice[2])):
        k = format(i, 'b').zfill(len(o_sublattice[2]))
        sum = 0
        for m in range(len(k)):
            sum += int(k[m])
        if sum == atom_num:
            pattern.append(k)
    return pattern


def unique(parent_sym_jun, pattern):
    """
    得られた究極の対称操作の辞書に基づいて、配列をユニークしていく

    parameters

    parent_sym_jun

    omomi4 は各空孔数ごとのcandidate集合

    retruns

    """
    pattern_copy = deepcopy(pattern)  # copy wo sakusei
    lis = set()

    for i in range(len(pattern)):  # through all candidate
        if pattern[i] in pattern_copy:  # kouho ga mada 生き残ってるかチェック
            for j in range(1, len(parent_sym_jun)):  # 置換操作について回す　
                d = dict()
                sta = ""
                for k in range(len(parent_sym_jun[0])):
                    d[parent_sym_jun[j][k]] = pattern[i][k]  # str辞書の作成
                # staの作成 sta is made from tikan[j]
                for r in range(len(parent_sym_jun[0])):
                    sta += d[r]
                if int(sta) is not int(pattern[i]):
                    if sta in pattern_copy:
                        pattern_copy.remove(sta)
                        lis.add(pattern[i])
    lis = list(lis)
    return lis


def sucell_unique(rutile, HNF, index, parent_sym):
    parent_lattice = get_superlattice(rutile, HNF, index)
    o_sublattice = get_o_superlattice(rutile, HNF, index)
    goalen = shin_get_permutation(o_sublattice, parent_sym, HNF)
    ds_list = list()
    for oxy_number in range(index*4):
        can = make_candidate(o_sublattice, oxy_number)
        auncko = unique(goalen, can)
        ds_list.extend(auncko)
    return ds_list


def unique_nov(rotational_permutations, mate_list):

    uniqued_mates = set()
    mate_dict = candidate_list_to_dict(mate_list)
    reduced_permutations = eliminate_identity_perm(rotational_permutations)

    for mate in mate_list:
        # 構造配列が生き残っているかどうか
        if mate_dict[mate]:
            # 置換操作について回す, 恒等操作を除く
            for perm in reduced_permutations:
                mate_prime = perm_on_mate(perm, mate)
                if int(mate) != int(mate_prime):
                    mate_dict[mate_prime] = False
                    uniqued_mates.add(mate)
    return list(uniqued_mates)


def perm_on_mate(perm, mate):
    # permutationをmate配列に作用させる
    sub_dict = {}
    ans = ''
    for site_before, site_after in perm.items():
        sub_dict[site_after] = mate[site_before]
    # answer mate wo create
    num_sites = len(sub_dict)
    for num in range(num_sites):
        ans += sub_dict[num]

    return ans


def eliminate_identity_perm(rotational_permutations):

    reduced_permutations = []
    for permutation in rotational_permutations:
        if is_identical_permuation(permutation):
            pass
        else:
            reduced_permutations.append(permutation)
    return reduced_permutations


def is_identical_permuation(permutaion):
    for i, j in permutaion.items():
        if i == j:
            continue
        return False
    return True


def candidate_list_to_dict(mate_list):
    mate_dict = {}
    for mate in mate_list:
        mate_dict[mate] = True
    return mate_dict
