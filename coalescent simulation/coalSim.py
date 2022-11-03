import numpy
import random
import pprint 
import scipy.stats as stats
#import lacData
#import alleleHist


"""
import A        # module named "A"

from A import B # module named "A", with member named "B"
                # member = individual variable or function 
import A.B      # module named "A.B"

import scipy.stats    #preferred way :)
scipy.stats.blah...

from scipy import stats #boo
stats.blah...

import scipy.stats as stats # also good
stats.blah...
"""

def next_coalescence(num_alleles, pop_size):
    """Returns a geometric random variable, the time to the next
    coalescence event (in generations)."""
    p = num_alleles * (num_alleles-1) / (4 * pop_size)
    return int(stats.geom.rvs(p)) 


def random_coalescent_tree(num_alleles, pop_size):
    """takes in population size and the number of alleles in the sample and 
    then creates a random, neutral coalescent tree with the following format:
    (node_num, left_tree, right_tree, tot_gens_before_present
    node_num and tot_gens_before_present are integers, left_tree and right_tree are tuples"""

    alleles_list = [(i,(),(),0) for i in range(num_alleles)]
    alleles_count = num_alleles 
    
    while len(alleles_list) >= 2:  
        #randomly choose 2 alleles for the first coalescence
        allele_one, allele_two = random.sample(alleles_list,2)

        #obtain time to coalescence 
        time = next_coalescence(len(alleles_list), pop_size)

        #update tree/build new node
        new_node = (alleles_count, allele_one, allele_two, time)
        alleles_count += 1 

        #update alleles_list
        alleles_list.append(new_node)
        alleles_list.remove(allele_one)
        alleles_list.remove(allele_two)

    return alleles_list[0]


pprint.pprint(random_coalescent_tree(5, 1000))




