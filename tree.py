from __future__ import print_function

def make(scores, comparator = min):
    sort_tup = lambda *x: tuple(sorted(x))
    def calculate_divergence_scores(distance_matrix):
        names = set(reduce(lambda x,y:x+y, distance_matrix, tuple()))
        return {name:sum([value for key, value in distance_matrix.iteritems() if name in key]) for name in names}
    def distance_to_own_node(taxa, distance_matrix):
        divergence_scores = calculate_divergence_scores(distance_matrix)
        names = divergence_scores.keys()
        lt, rt = sort_tup(*taxa)
        lt_dist = 0.5*distance_matrix[(lt,rt)] + (divergence_scores[lt] - divergence_scores[rt]) / (2*(len(names) - 2))
        rt_dist =  0.5*distance_matrix[(lt,rt)] + (divergence_scores[rt] - divergence_scores[lt]) / (2*(len(names) - 2))
        return {lt:lt_dist, rt:rt_dist}
    def distance_to_node(tax, node, distance_matrix):
        return (distance_matrix[sort_tup(tax, node[0])] + distance_matrix[sort_tup(tax, node[1])] - distance_matrix[node]) / 2 
    def modified_distance(lr, distance_matrix):
        divergence_scores = calculate_divergence_scores(distance_matrix)
        names = divergence_scores.keys()
        l,r = sort_tup(*lr) 
        return distance_matrix[lr] - (divergence_scores[l] + divergence_scores[r])/(len(names) -  2)

    scores = dict(scores)
    nodes = []
    distance_matrix = {sort_tup(*key): scores[key] for key in scores} 
    while len(distance_matrix) > 1:
        q_matrix = {lr:modified_distance(lr, distance_matrix) for lr in distance_matrix} #calc modified distance for each pair
        q_matrix_names = set(reduce(lambda x,y: x+y, q_matrix, tuple())) #all the nodenames
        new_node = comparator(q_matrix, key = lambda x:q_matrix[x]) #select the most appropriate new node
        nodes.append((new_node,distance_to_own_node(new_node, distance_matrix))) #add the new node to the results
        for key in q_matrix_names:
            if key in new_node:
                continue
            distance_matrix[sort_tup(key, new_node)] = distance_to_node(key, new_node, distance_matrix)
            for child in new_node: 
                del q_matrix[sort_tup(key, child)]
                del distance_matrix[sort_tup(key, child)]
        del distance_matrix[new_node]
    new_node = list(distance_matrix)[0]
    nodes.append((new_node, distance_matrix[new_node]))
    return nodes


def to_newick(tree):
    def replace(cur_node, branch_lengths):
        if isinstance(branch_lengths, dict):
            results = list(cur_node) + [branch_lengths[child_node] for child_node in cur_node]
            #this gives us a list formatted like [A,5,B,10] or [A,5, (B,C),7] or [(A,B),5,[C,D),7]
            outstr = '({0}:{2:.9f}, {1}:{3:.9f})'
        else: #it is the top-node
            results = list(cur_node + (branch_lengths,))
            outstr = '({0},{1}:{2:.9f})'
        if str(cur_node).count('(') != 1: #this has 1 or more non-leaf child-nodes
            for pos, child_node in enumerate(cur_node): #
                if str(child_node).count('(') == 0: #this is a leaf-node
                    continue
                for node ,node_branch_lengths in tree: 
                    if child_node == node: #this is the node to which the branch points
                        results[pos] = str(replace(node, node_branch_lengths))
                        break
        return outstr.format(*results)
    tree = tree[:]
    top, topdist = tree[-1]
    return replace(top, topdist)

def print_tree(tree):
    for node, branch_lengths in sorted(tree, key = lambda x:str(x[0]).count('(')):
        if hasattr(branch_lengths, '__iter__'):
            print(node)
            for child_node in branch_lengths:
                print(child_node, branch_lengths[child_node])
        else:
            print(branch_lengths)
        print('----------------')
