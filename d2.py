#----default libraries ---- #

from __future__ import division, print_function
from __main__ import eprint
import math
import itertools
from collections import namedtuple, Counter

__as_aminoacids__ = False
__with_substitutes__ = False

substitions_letters = set('RYSWKMBDHVN')
substitution_table = {
    'A':'A',
    'C':'C',
    'G':'G',
    'T':'T',
    'R':'AG',
    'Y':'CT', 
    'S':'GC', 
    'W':'AT',
    'K':'GT',
    'M':'AC',
    'B':'CGT',
    'D':'AGT',
    'H':'ACT',
    'V':'ACG',
    'N':'ACGT'
                        
}
flipsub_table = {
    'A':'ARWMDHVN',
    'C':'CYSMBHVN',
    'G':'GRSKBDVN',
    'T':'TYWKBDHN',
}

subword_cost = {}



class SeqHolder(namedtuple('SeqHolder', 'adj_len letter_dist word_count')):
    '''adj_len is the adjusted sequence_length,adj_lenght = sequence_length - word_length'''
    def word_prob(self, word):
        product = lambda iterable: reduce(lambda x,y:x*y, iterable)
        return product((self.letter_dist[letter] for letter in word))
    def normalized_frequency(self, word):
        return float(self.word_count[word]) - self.adj_len * self.word_prob(word)
    def log_normalized_frequency(self, word):
        #this is another way to have a zero-centered, normally distributed frequency function, which is required for shepp distance
        if (not self.count(word)): #not observed
            if (not self.word_prob(word)): #not expected either
                return 0 #the number of expected matched the number of observed, so return 0
            else: 
                return math.log(float(self.count(word) + 1) / (self.adj_len * self.word_prob(word)), 2) 
                #cannot take a log of 0, so instead we take a log of the minimum observable
        else:
            ratio = float(self.count(word)) / (self.adj_len * self.word_prob(word))
            return math.log(ratio, 2)
    def setup_substcount(self, possible_words):
        words_copy = self.word_count.copy()
        for word in possible_words:
            self.word_count[word] = sum((float(words_copy[subst_word]) / subst_value
                                for subst_word, subst_value in subword_cost[word]))
    def count(self, word):
        #if __with_substitutes__:
        #    return sum((float(self.word_count[subst_word]) / subst_value
         #                       for  subst_word, subst_value in create_substitutes(word)))
        return self.word_count[word]

def create_substitutes(word):
    product = lambda iterable: reduce(lambda x,y: x*y, iterable, 1)
    subst_words = map(''.join, itertools.product(*map(lambda letter: flipsub_table[letter], word)))
    subst_values = map(lambda word: product((len(substitution_table[letter]) for letter in word)), subst_words) 
    return zip(subst_words, subst_values)

def setup_substitutes(possible_words):
    global subword_cost #strictly not needed, but for clarities sake
    product = lambda iterable: reduce(lambda x,y: x*y, iterable, 1)
    for word in possible_words:
    #what matches do we have
        subst_words = map(''.join, itertools.product(*map(lambda letter: flipsub_table[letter], word)))
    #how many words does this match represent
        subst_values = map(lambda word: product((len(substitution_table[letter]) for letter in word)), subst_words) 
        subword_cost[word] = zip(subst_words, subst_values)

def setup_holder(seq, word_length, alphabet):
    get_words = lambda seq, length: [seq[x:x+length] for x in range(len(seq) - length + 1)]
    seq_len = float(len(seq))
    adj_len = seq_len - word_length
    letter_dist = {letter:seq.count(letter)/seq_len for letter in (substitution_table.keys() if __with_substitutes__ else
                                                                    alphabet)}
    #print letter_dist
    #raw_input()
    word_count = Counter(get_words(seq, word_length))
    return SeqHolder(adj_len, letter_dist, word_count)


def d2_distance(score, seqholders, possible_words):
    product = lambda iterable: reduce(lambda x,y:x*y, iterable)
    pyth = lambda iterable:math.sqrt(reduce(lambda x,y:x + y**2.0, iterable, 0)) #pythogarean distance,    sub = 
    pyths = [pyth((seq.count(word) for word in possible_words)) for seq in seqholders]
    sub = product(pyths)
    
    try:
        return 0.5 * (1 - float(score)/sub)
    except ZeroDivisionError:
        print("fail")
        for pos, seq in enumerate(seqholders):
            print(seq, "all zeros {}".format(not any(seq.word_count.values())))
        print (pyths)
        print(create_substitutes("AAAA"))
        raise

def dshep_distance(score, seqholders, possible_words, normalize = SeqHolder.normalized_frequency):
    '''d2s as in song et al pg 347
    ~x[w] = normalized_frequency count of word w in sequence x 
    pyth = pythogarean distance  = (a^2+b^2)^0.5
    ds2 =  0.5 * ( 1 - (DS2 
                        / 
        ((~x[w]^2)/pyth(~x[w],~y[w]) for w in words) *  ((~y[w]^2)/pyth(~x[w],~y[w]) for w in words)
        ))
    
    '''
    product = lambda iterable: reduce(lambda x,y:x*y, iterable)
    pyth = lambda iterable:math.sqrt(reduce(lambda x,y:x + y**2.0, iterable, 0)) #pythogarean distance,
    pyth_count = lambda seqholders, word: pyth((normalize(seq, word) for seq in seqholders)) 
    #pythogorean distance of normalized counts of given word in given counters
    return 0.5 * (1 - float(score)                                                  \
            /                                                                       \
        product((                                                                   \
            sum(                                                                    \
                (normalize(seq, word)**2.0                         \
                                /                                                   \
                pyth_count(seqholders, word) for word in possible_words)              \
                ) ** 0.5                                                            \
            for seq in seqholders)))


def dstar_distance(score, seqholders, possible_words):
    product = lambda iterable: reduce(lambda x,y:x*y, iterable)
    correction = product((                                                                  \
        sum(                                                                    \
            (seq.normalized_frequency(word) ** 2.0                       \
                /                                                               \
            (seq.adj_len * seq.word_prob(word)) for word in possible_words)      \
            ) ** 0.5                                                            \
        for seq in seqholders))
    return 0.5 * (1 - float(score)/correction)

def dstar_distance_debug(score, seqholders, possible_words):
    product = lambda iterable: reduce(lambda x,y:x*y, iterable)
    def cor_func(num, seq):
        for word in possible_words:
            try: 
                yield seq.normalized_frequency(word) ** 2.0 \
                    / seq.adj_len * seq.word_prob(word) 
            except ZeroDivisionError:
                print(num ,seq.adj_len, seq.word_prob(word), word)
                raise
        
    correction = product((                                                                  \
        sum(cor_func(num, seq)) ** 0.5                                                            \
        for num, seq in enumerate(seqholders)))
    return 0.5 * (1 - float(score)/correction)


def score_pair(seqholders, mode = 'D2', possible_words = None, as_distance = False):
    product = lambda iterable: reduce(lambda x,y:x*y, iterable)
    def d2_score(seqholders, word):
        return product((seq.count(word) for seq in seqholders))
    
    def dstar_score(seqholders, word):
        return product((seq.normalized_frequency(word) for seq in seqholders))  \
                /                                                       \
            math.sqrt(
                product((seq.word_prob(word) * seq.adj_len for seq in seqholders))
                )
    def dshep_score(seqholders, word, normalize = SeqHolder.normalized_frequency):
        return \
            product((normalize(seq, word) for seq in seqholders))  \
                /                                                               \
            math.sqrt(
                sum((normalize(seq, word)**2 for seq in seqholders))
                    )
        
    stat =  {
         'D2' : d2_score,
         'D*2': dstar_score,
         'DS2': dshep_score,
         'DS2L':lambda x,y : dshep_score(x,y, normalize = SeqHolder.log_normalized_frequency)
            } 
    assert mode in stat
    assert possible_words is not None
    
    
    score = sum((stat[mode](seqholders, word) for word in possible_words))
    
    if as_distance:
        return {
            'D*2':dstar_distance,
            'DS2':dshep_distance,
            'DS2L':lambda x,y,z : dshep_distance(x,y,z, normalize = SeqHolder.log_normalized_frequency),
            'D2':d2_distance,
            }[mode](score, seqholders, possible_words)
    else:
        return score    
    

def make_alphabet(sequences):
    '''infer the alphabet from the supplied sequences'''
    global __with_substitutes__
    allchars = set(''.join(sequences))
    if (not __as_aminoacids__)  and substitions_letters & allchars: #we have a sequence with substitions
        __with_substitutes__ = True
        return ''.join(allchars - substitions_letters)
    else:
        return ''.join(allchars) 
    
def make_possible_words(alphabet, word_length):
    '''returns list of all possible word_length words for a given alphabet'''
    return [''.join(letters) for letters in itertools.product(alphabet, repeat = word_length)]


def word_count_sequences(sequences,word_length,alphabet = None):
    '''returns seqholder objects for given sequences and k-mer length
    if no alphabet is supplied, it will infer one'''
    if alphabet is None: 
        alphabet = make_alphabet(sequences)
    return [setup_holder(seq,word_length, alphabet) for seq in sequences]

def align(named_sequences, word_length,  mode = 'D2', as_distance = False, pair_size = 2):
    '''takes a dict of {seqname:seq} or list of (seqname, seq) and 
    returns list of tuples of the score between each pair as 
    ((seqname_a, seqname_b), score)
    
    mode:
        D2, D*2 for D-Star, DS2 for D-shepherd, DS2L for log D-shepher
    as_distance:
        calculates ds2 for DS2 or d*2 for D*2
    pair_size:
        untested, 
    
    '''
    if isinstance(named_sequences, dict):
        names, sequences = zip(*named_sequences.iteritems())
    elif isinstance(named_sequences, list):
        names, sequences = zip(*named_sequences)
    else:
        raise RuntimeError("named_sequences must be a list or dict")
    alphabet = make_alphabet(sequences)
    seqholders = word_count_sequences(sequences, word_length, alphabet)
    possible_words = make_possible_words(alphabet, word_length)
    
    if __with_substitutes__: #we setup the hashtable for substitutes to improve performance here, rather than recalculating every time
        setup_substitutes(possible_words)
        for seqholder in seqholders:
            seqholder.setup_substcount(possible_words) #adjust the wordcount according to the hastable created above
        
    aligns = []
    for pair_numbers in itertools.combinations(range(len(sequences)), pair_size):
        pair_names = tuple([names[nr] for nr in pair_numbers])
        groupholders = tuple([seqholders[nr] for nr in pair_numbers])
        alignscore = score_pair(groupholders, mode = mode,
                                possible_words = possible_words,
                                as_distance = as_distance,
                                )
        aligns.append((pair_names, alignscore))
    return aligns





















