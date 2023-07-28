from Bio import SeqIO
import re
import os
import sys

SEQUENCES_PATH = sys.argv[1] ##
##SEQUENCES_PATH = 'FNA/'

if SEQUENCES_PATH.endswith("/"):
    SequencesPath = re.sub('/', '', SEQUENCES_PATH)
    SequencesDir = str(SEQUENCES_PATH)
else:
    SequencesPath = SEQUENCES_PATH
    SequencesDir = str("".join ([SEQUENCES_PATH,'/']))
    
Orthologues = [x for x in os.listdir(SequencesDir) if x.endswith(".fasta")]

for Sequence in Orthologues:
    FNA = str("".join ([SequencesDir,Sequence]))
    FNAOut = str("".join ([SequencesDir,Sequence,".phy"]))
    records = SeqIO.parse(FNA, "fasta")
    count = SeqIO.write(records, FNAOut, "phylip")
