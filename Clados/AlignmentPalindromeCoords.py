from Bio import SeqIO
import re
import os
import sys

SEQUENCES_PATH = sys.argv[1]
#SEQUENCES_PATH = 'gbff/Only_ORTHOLOGUES/' ## Path de la carpeta con los ortólogos.
pattern = sys.argv[2]
pattern = re.sub('G', '[gG][-]*', pattern)
pattern = re.sub('C', '[cC][-]*', pattern)
pattern = re.sub('T', '[tT][-]*', pattern)
pattern = re.sub('A', '[aA][-]*', pattern)
pattern = pattern[:-4]
SPP = sys.argv[3]
output_file = 'Orthologues_Palindrome_sites.txt' ## Nombre del archivo de salida
output = open (output_file, 'w') ## Abrimos el archivo de salida
output.write('FILE\tPAL\tSTART\tEND\n')

if SEQUENCES_PATH.endswith("/"): ## Esta parte la pongo por si corro este codigo en la terminal
    SequencesPath = re.sub('/', '', SEQUENCES_PATH) ## De este modo evito errores en el argumento path 
    SequencesDir = str(SEQUENCES_PATH)
else:
    SequencesPath = SEQUENCES_PATH
    SequencesDir = str("".join ([SEQUENCES_PATH,'/']))
Orthologues = [x for x in os.listdir(SequencesDir) if x.endswith(".phy")] ## creo un arreglo con todos los ortólogis de la carpeta
#pattern = '[tT][-]*[cC][-]*[gG][-]*[aA][-]*[tT][-]*[cC][-]*[gG][-]*[aA]'
j=0
k=0
for Orthologue in Orthologues:
    FNA = str("".join ([SequencesDir,Orthologue])) ## creo el path completo con el path de la carpeta de ortólogos y el path del ortólogo
    for record in SeqIO.parse(open(FNA),'phylip'):
        Seq = str(record.seq)
        Spp = record.description
        i = 0
        if Spp == SPP:#PCC_6303_PCC_6303 336-3_336-3
            for match in re.finditer(pattern, Seq):
                i += 1
                j += 1
                site = match.group()
                start = match.span()[0]
                end = match.span()[1]
                if len(site)==8:
                    k += 1
                    print ('{}\t{}\t{}\t{}:{}'.format(i,Orthologue,site,start+1,end))
                    output.write('{}\t{}\t{}\t{}\n'.format(Orthologue,site,start+1,end))
    print ('_________________________________________________________________________\n')
output.close()
print ("TOTAL: {} sitios".format(k))
