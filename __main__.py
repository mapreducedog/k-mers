#---default libraries----#
from __future__ import print_function
import sys
#declare here so other d2 and tree can use this
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


#-----own-libraries-----#
import d2
import tree
import fasta_paml_io as fpio


helpstr =  ''' d2.py, a d2 matching and neighbour joining script by MDP.
        usage: __main__.py [options] path/to/file 
        options
        -mode ['D2', 'D*2','DS2', DS2L']
            matching algorithm to use, normal D2, 
            star or shepherd, default D*2
        -length [1,2...10]
            length of k-mers to compare, default 4}
            
        -NJ [1/0]
            perform a neighbour join, default 1
            
    '''


def parse_args(args):
    '''
    takes list like ['-mode', 'D*2', -NJ,1], list can be empty,    
    default args:
    -mode:'D*2'
    -NJ:1
    -length:4
    -debug:0
    ''' 
    
    validopts = {'-mode':['D2', 'D*2', 'DS2', 'DS2L'], '-NJ':xrange(2),'-length':xrange(1,100), '-debug':xrange(2)}
    opts = {'-mode':'D*2', '-NJ':1, '-length':4, '-debug':0}
    if args:
        givenopts = dict(zip(args[0::2],args[1::2]))
        for key,value in givenopts.iteritems():
            if key in {'-NJ', '-length', '-debug'}:
                try:
                    value = int(value)
                except ValueError:
                    raise KeyError("Invalid input: {} is not a valid argument for {}".format(value, key))
            if key not in validopts:
                raise KeyError("Invalid input: {} is not an option".format(key))
            if value not in validopts[key]:
                raise KeyError("Invalid input: {} is not a valid argument for {}".format(value, key))
            opts[key] = value
    return opts


def perform_align(filename, opts):
    sqs = fpio.get_sequences(filename)
    align = d2.align(sqs, opts['-length'], mode = opts['-mode'], as_distance = opts['-NJ'])
    eprint("passed align!")
    #align = align_groups(sqs, lambda x:d2_score(opts['-length'], x, opts['-mode'], as_distance = True))
    # if opts['-debug']:
        # return align
    # dist_align = pair_dist(align)
    if opts['-NJ']:
        mytree = tree.make(align, min)
        return tree.to_newick(mytree)
    else:
        return sorted(align, key = lambda x:x[-1])

def main():
    if {'-help','--help','-h','--h'} & set(sys.argv) or len(sys.argv) == 1:
        return helpstr
    filename = sys.argv[-1]
    #if len(sys.argv) > 2:
    args = sys.argv[1:-1]
    #else:
    #    args = []
    opts = parse_args(args)
    return perform_align(filename, opts)
if __name__ == '__main__':
    print(main())

