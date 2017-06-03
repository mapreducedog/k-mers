def fasta_parser(filename, format):
    sequences = {}
    name = ''
    with open(filename) as sequencefile :
        for line in sequencefile.readlines():
            line = line.strip()
            if line.startswith('>'):
                name = line[1:]
                sequences[name]  = ''
            elif line: #its not empty
                sequences[name] += line.strip()
    if format is dict:
        return sequences
    else:
        return list(sequences.itervalues())
        
def paml_parser(filename, format):
    with open(filename) as sequencefile: 
        if format is dict:
            return {line.partition(' ')[0].strip() : ''.join(line.strip().split()[1:]) for line in sequencefile.readlines() if line[0].isalpha()}
        else:
            return [''.join(line.strip().split()[1:]) for line in sequencefile.readlines() if line[0].isalpha()]
        

def get_sequences(filename,  format = dict):
    '''reads a .paml or .fasta 
    and  returns either a dict with {seq_name:seq} or list of (seq_name,seq)
    raises RunTimeError if not a .paml or .fasta'''
    if filename.endswith('.fasta'):
        return fasta_parser(filename, format)
    if filename.endswith('.paml'):
        return paml_parser(filename, format)
    raise RuntimeError("Couldnt parse file, unknown extension {}, expected .fasta or .paml".format(filename.split('.')[-1]))

