'''
Generate the biased ab initio forward folding input files

The command will be:
    
python lowrms_frags_topN.py -frag_qual frag_qual3.dat -ntop 3 
    -fullmer input.200.3mers -out input.3t200.3mers

'''
from argparse import ArgumentParser

# Add user interface arguments:
parser = ArgumentParser(description = 'Get the lowest rmsd fragments among top N')
parser.add_argument('-frag_qual', type = str, help = 'fragment quality file storing RMSD and score ranking')
parser.add_argument('-fullmer', type = str, help = 'fragment file')
parser.add_argument('-out', type = str, help = 'output fragment file')
parser.add_argument('-ntop', type = int, default = 3, help = 'How many fragments of lowest RMSD to take')
args = parser.parse_args()

frag_qual_file = args.frag_qual
fullmer = args.fullmer
outfile = args.out
ntop = args.ntop

# Read in the frag.fsc file to store all identifying info + rmsd
filein = open(frag_qual_file)
unranked_frags = {}
nfrag = 0

for line in filein:
    if line.startswith('#'):
        nfrag = 0
        continue
    
    else:
        nfrag += 1
        # Get position of fragment on input pdb
        pos = int(line.split()[0])
        
        # Get rmsd of fragment
        rmsd = float(line.split()[-3])
        
        if nfrag == 1:
            unranked_frags[pos] = {}
        
        unranked_frags[pos][nfrag] = rmsd
        
filein.close()

# Find the N lowest RMSD fragments:
ranked_frags = {}

for pos in unranked_frags.keys():
    frag_list = list(unranked_frags[pos].keys())
    frag_list = sorted(frag_list, key=lambda x: unranked_frags[pos][x])
    top_n_list = frag_list[:ntop]
    ranked_frags[pos] = top_n_list

# Take the N best fragments from frag_id_num.200.9mers
filein = open(fullmer)
fileout = open(outfile, 'w')

for line in filein:
    if 'position' in line:
        pos = int(line.split()[1])
        if pos > 1:
            fileout.write('\n')
        newline = line.replace('200','%d'%ntop)
        fileout.write(newline)
        prev_frag = 0
        count = 0
        nfrag = 0
            
    elif line == '\n':
        nfrag += 1
            
    elif len(line.split()) > 1:
        if nfrag in ranked_frags[pos]:
            this_frag = nfrag
            if this_frag != prev_frag:
                count += 1
                fileout.write('\n')
                prev_frag = this_frag
            fileout.write(line)