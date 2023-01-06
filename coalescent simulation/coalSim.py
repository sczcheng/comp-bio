import numpy
import random
import pprint
import scipy.stats as stats


def next_coalescence(num_alleles, pop_size):
    """returns a geometric random variable, the time to the next
    coalescence event (in generations)."""
    p = num_alleles * (num_alleles - 1) / (4 * pop_size)
    return int(stats.geom.rvs(p))


def random_coalescent_tree(num_alleles: int, pop_size: int):
    """takes in population size and the number of alleles in the sample and
    then creates a random, neutral coalescent tree with the following format:
    (node_num, left_tree, right_tree, tot_gens_before_present)
    node_num and tot_gens_before_present are integers, left_tree and right_tree are tuples"""

    alleles_dict = {i: (i, (), (), 0) for i in range(num_alleles)}
    alleles_count = num_alleles

    while len(alleles_dict) >= 2:
        # choose 2 alleles to coalesce
        (label_1, subtree_1), (label_2, subtree_2) = random.sample(
            list(alleles_dict.items()), 2
        )

        # obtain time to coalescence
        time = next_coalescence(len(alleles_dict), pop_size)

        # update tree and alleles_dict
        alleles_dict[alleles_count] = (alleles_count, subtree_1, subtree_2, time)
        alleles_count += 1
        del alleles_dict[label_1]
        del alleles_dict[label_2]

    return next(iter(alleles_dict.values()))


##### testing #####
# pprint.pprint(random_coalescent_tree(5, 500))


def convert_to_gens_per_branch(tree, parent_gens_before_present):
    """takes in a tree with where each node has coalescence in units of generations,
    and then returns a tree in the following format (node_num, left_tree, right_tree, branch_gens)
    where branch_gens is the number of generations on the branch leading on node_num. The input parameter
    parents_gens_before_present is the generation time of the parent coalescence at the start of branch_gens
    leading to node_num
    """

    if tree == ():
        return tree

    return (
        tree[0],
        convert_to_gens_per_branch(tree[1], tree[-1]),
        convert_to_gens_per_branch(tree[2], tree[-1]),
        parent_gens_before_present - tree[-1],
    )


##### testing #####
# tree_1 = (4, (2, (), (), 0), (3, (1, (), (), 0), (0, (), (), 0), 20), 54)
# time_of_final_coalescence=tree_1[-1]
# pprint.pprint(convert_to_gens_per_branch(tree_1,time_of_final_coalescence))


def get_branch_proportion_list(tree):
    """returns a list containing branch numbers where the number of times each branch number appears
    is the number of generations on that branch"""

    if tree == ():
        return []

    return (
        [tree[0]] * tree[-1]
        + get_branch_proportion_list(tree[1])
        + get_branch_proportion_list(tree[2])
    )


##### testing #####
# tree=(4, (2, (), (), 4), (3, (1, (), (), 1), (0, (), (), 1), 2), 5)
# pprint.pprint(get_branch_proportion_list(tree))


def assign_muts_to_branch(tree, num_muts):
    """takes in a tree and the number of mutations to lay down on the tree.
    Then, randomly assigns numbered mutations to the branches and returns a list of lists. The index in the
    outer list correspond to branch number, whereas the inner list at that index contains the mutations that occur on a branch
    """

    branch_muts_proportion = get_branch_proportion_list(tree)

    # tree[0]+1 is also the total number of nodes
    assigned_mutations = [[] for _ in range(tree[0] + 1)]

    for mut_id in range(num_muts):
        node_id = random.choice(branch_muts_proportion)
        assigned_mutations[node_id].append(mut_id)

    return assigned_mutations


##### testing #####
# tree=(4, (2, (), (), 4), (3, (1, (), (), 1), (0, (), (), 1), 2), 0)
# pprint.pprint(assign_muts_to_branch(tree, 4))


def create_seqs(tree, branch_mut_list, seq):
    """this function recursively traverses a tree applying mutations on branches. The input branch_mut_list is the output of the function assign_muts_to_branch.
    The input seq is the starting sequence inherited at the root of the tree. The function returns a list of sequences, which are
    tuples containing mutations found on them.
    """
    if not tree:
        return []

    node_num, left_tree, right_tree, _ = tree
    mut = (*seq, *branch_mut_list[node_num])

    left_n_right = [
        *create_seqs(left_tree, branch_mut_list, mut),
        *create_seqs(right_tree, branch_mut_list, mut),
    ]

    if not branch_mut_list[node_num]:
        return left_n_right
    else:
        return [mut, *left_n_right]


##### testing #####
# tree = (4, (2, (), (), 4), (3, (1, (), (), 1), (0, (), (), 1), 2), 0)
# branch_mut_list = [[3], [], [0, 2], [1], []]
# pprint.pprint(create_seqs(tree, branch_mut_list, ()))


def coal_sim(num_alleles, pop_size, num_muts):
    """takes in number of alleles, population size, and number of mutations and performs a coalescent simulation
    This function outputs a list of sequences (haplotypes) created when mutations are put on a random neutral tree"""

    rct = random_coalescent_tree(num_alleles, pop_size)
    tree = convert_to_gens_per_branch(rct, rct[-1])

    return create_seqs(tree, assign_muts_to_branch(tree, num_muts), ())


##### testing #####
# pprint.pprint(coal_sim(5, 200, 8))
