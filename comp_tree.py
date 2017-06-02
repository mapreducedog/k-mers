from __future__ import print_function
import re
import sys
import glob
__doc__ = '''comp_tree.py
compares newick for the number of shared nodes
usage : comp_tree [<glob_for_files>]
if no path is supplied, inspect the current working directory for all ".txt" files, 
otherwise globs as supplied argument
for example 
python comp_tree ~/*.txt
or
python comp_tree ~/*.tree
'''


#example = r"(SEQCCC,(SEQDDD:0.060489303, ((SEQBBB:0.118115857, (SEQAAA:0.104426401, (SEQGGG:0.150950683, SEQHHH:0.137104714):0.050034382):0.011685270):0.007855961, (SEQEEE:0.045948759, SEQFFF:0.044118186):0.032791438):0.006207807):0.074406265)"
example = r"(((A,B),(C,D)), E)"
match_values = re.compile(r"(:[\d\.]+)|;") #replaces values and ";"
match_names = re.compile(r"([a-zA-Z0-9]+)")
match_lbrace = re.compile(r"\(")
match_rbrace = re.compile(r"\)")
matches_re = [match_values, match_names, match_lbrace, match_rbrace]
substitutes = [ '', r'"\1"', "[", "]"]
methods = ["D2", "DS2", "DS2L", "Answer"]
def compare_nodes(query_nodes, answer_nodes):
    #print(query_nodes,"\n", answer_nodes, "\n", query_nodes & answer_nodes,"\n\n")
    return len(query_nodes & answer_nodes)


def substitute(instr, match_sub):
    match, sub = match_sub
    return match.sub(sub, instr)

def normalize_tree(instr_tree):
    str_tree = reduce(substitute, zip(matches_re, substitutes), instr_tree)
    return eval(str_tree)
def __get_nodes__(ltree):
    '''create a flat list of nodes,
    #that is to say ((((A,B), C), (D,E)),F ) or
    (
        (
            (
                (A,B), 
            C), 
        (D,E)
        )
    ,F )
    becomes  : (A,B),(A,B,C),(D,E), (A,B,C,D,E), (A,B,C,D,E,F)
    '''
    if isinstance(ltree, str): #this is a leaf
        return (frozenset(), frozenset([ltree,]))
    total, latest = reduce(
        #(total, latest), (childtotal, childlatest) -> (total | childtotal, latest | childlatest)
        lambda cur_tot_lat, child_tot_lat : (cur_tot_lat[0] | child_tot_lat[0], cur_tot_lat[1] | child_tot_lat[1]),
        map(__get_nodes__, ltree),
        (frozenset(), frozenset()))
    total |= {latest, } #add its own group to the total
    return total, latest
def get_nodes(ltree):
    total, top = __get_nodes__(ltree)
    return total - {top,} #remove the headnode from the list, everyone will share this
def preprocess_file(filename):
    with open(filename) as infile:
        lines = map(lambda x: x.strip(), filter(lambda x:x.startswith('('), infile))
    return lines
def rank_trees(treeslist):
    '''expects a list of trees of which the last tree is the "answer_tree, returns a list of total matching nodes between those trees and the answer trees, (the last column is thus the total number of (non-terminal, non-head) nodes in the answer_tree'''
    answer_tree = treeslist[-1]
    answer_nodes = get_nodes(answer_tree)
    return map(lambda query_tree : compare_nodes(get_nodes(query_tree), answer_nodes), treeslist)
    

def handle_files():
    scores = []
    trees_ensemble = []
    globstr = sys.argv[1] if len(sys.argv) > 1 else '*.txt'
    for filename in sorted(glob.glob(globstr)):
        trees_ensemble.append(map(normalize_tree, preprocess_file(filename)))
    scores = map(rank_trees, trees_ensemble)
    return scores


if __name__ == "__main__":
    if '-h' in sys.argv:
        print(__doc__)
    else:
        print(*methods)
        scores = handle_files()
        map(print, scores)
    
