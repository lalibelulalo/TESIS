from Bio import SeqIO
import re
import os
import sys
import shutil

SEQUENCES_PATH = sys.argv[1] #'Orthologues_GCGATCGC_336_3/' ## Path de la carpeta con los ortólogos.
TaxonNumber = int(sys.argv[2])
#output_file = 'Orthologues_Palindrome_sites.txt' ## Nombre del archivo de salida
#output = open (output_file, 'w') ## Abrimos el archivo de salida

if SEQUENCES_PATH.endswith("/"): ## Esta parte la pongo por si corro este codigo en la terminal
    SequencesPath = re.sub('/', '', SEQUENCES_PATH) ## De este modo evito errores en el argumento path 
    SequencesDir = str(SEQUENCES_PATH)
else:
    SequencesPath = SEQUENCES_PATH
    SequencesDir = str("".join ([SEQUENCES_PATH,'/']))
Orthologues = [x for x in os.listdir(SequencesDir) if x.endswith(".fna")] ## creo un arreglo con todos los ortólogis de la carpeta

for Orthologue in Orthologues:
    i = 0
    FNA = str("".join ([SequencesDir,Orthologue])) ## creo el path completo con el path de la carpeta de ortólogos y el path del ortólogo
    FAA = re.sub('.fna', '.faa', FNA)
    OrthologueA = re.sub('.fna', '.faa', Orthologue)
    for record in SeqIO.parse(open(FNA),'fasta'):
        i += 1
    #print ('{}\t{}'.format(Orthologue,i))
    if i == TaxonNumber:
        shutil.copyfile(FNA,str("".join (['Only_ORTHOLOGUES/',Orthologue]))) ## si hasy unicamente 6 entradas entonces no hay paralogos
        shutil.copyfile(FAA,str("".join (['Only_ORTHOLOGUES/',OrthologueA])))
    else:
        shutil.copyfile(FNA,str("".join (['PARALOGUES/',Orthologue]))) ## Si hay mas de 6 entradas hay parálogos
        shutil.copyfile(FAA,str("".join (['PARALOGUES/',OrthologueA])))
