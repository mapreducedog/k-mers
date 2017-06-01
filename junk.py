#--default libraries---#
import itertools

def depth(branch):
    '''determines the max_depth of an multidimensional list'''
    if hasattr(branch, '__iter__') and not isinstance(branch, str):
        return max([depth(sub)+1 for sub in branch])
    else:
        return 0

def normalize(scores, inverted = False):
    numbers = [x[-1] for x in scores]
    top,bottom = max(numbers), min(numbers)
    if inverted:
        scale = lambda x: 1 - (float(x) - bottom)/(top - bottom)  
    else:
        scale = lambda x: (float(x) - bottom)/(top - bottom) 
    return [x[:-1]+(scale(x[-1]),) for x in scores]

        
def starting_points(sequence):
    return {sequence[:x] for x in range(1, len(sequence) + 1)}

def refmap(sequence):
    {sequence[x:x+y]:sequence.find(sequence[x:x+y]) for x in range(len(sequence)) for y in range(1, len(sequence[x:]) + 1)}
def unique(sequence, subseq_start, subseq_length):
    pos = 0
    subpos = 0
    subseq = sequence[subseq_start:subseq_start+subseq_length]
    while pos <=len(sequence) - subseq_length: #if the startign position is larger than the length, you can stop
        for subpos, char in enumerate(sequence[pos:pos+subseq_length]):
            if subseq[subpos] != char: #mismatch, skip to the 
                pos += subpos + 1
                break
        else: #It was a complete match
            return False
        if subseq_start < pos < subseq_start+subseq_length: #In the original sequence region
            pos = subseq_start+subseq_length
    return True
    
def find_shustrings(sequence):
    minlength = len(sequence) - 1
    shustrings = set()
    for start, letter in enumerate(sequence):
        length = 1
        while length <= minlength:
            subseq = sequence[start:start+length]
            if sequence.find(subseq) == sequence.rfind(subseq): #this is unique
                if length < minlength: #this is shorter than all previous unique subseq
                    minlength = length 
                    shustring = {subseq} #dispose of previous
                else:
                    shustring.add(subseq) #add to previous
            length += 1
    return shustrings


def align_groups(named_sequences,  function, groupsize = 2,  sort = None, match_self = False):
    '''takes a dic of {seq_name:seq} and performs the supplied function on each combination of sequences,
        returns a list of tuples of ((seqname_1, seqname_2), result) 
    options:
        groupsize:
            whether to perform function on duplets, triplets etc,
        sort:
            what order the return list should be sorted in
        match_self:
            whether function should be applied on (A,A) as well as (A,B) (combination with replacement)
    '''
    
    if sort is None:
        sort = {'reverse':True,  'key':lambda x: x[-1]} #sort from high scores to low scores
    else: 
        sort = sort.copy() #accounting for mutability
        if 'key' in sort: #make the sort function look at the score, rather than name
            key = sort['key'] 
            sort['key'] = lambda x: key(x[-1]) 
        else:
            key['key'] = lambda x: x[-1]
         

    if match_self:  #the comparison function also compares self with self.
        make_pairs = itertools.combinations_with_replacement
    else:
        make_pairs = itertools.combinations

#                                             |  sorts into ((seqname_A, seqname_B), (seq_A, seq_B) for each pair 
#                                             |    as it receives ((seqname_A, seq_A), (seqname_B, seq_B)) from make_pairs
#                                             v                                                   
    pair_names_seqs = itertools.imap(lambda pair: zip(*pair),
                        make_pairs(named_sequences.iteritems(), groupsize))
    
    return sorted(\
            map(lambda (names, seqs): (names, function(seqs)), pair_names_seqs),
            **sort)
        #return sorted([(group,  function([named_sequences[x] for x in group]))  for group in itertools.combinations(named_sequences.keys(),  groupsize)],  **sort)


