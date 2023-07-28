from Bio import SeqIO
import re
import os
import sys

SEQUENCES_PATH = sys.argv[1] #'RECONSTRUCCIONES/' #'Only_ORTHOLOGUES/' #
pattern = sys.argv[2] #'GCGATCGC' #
SPP = sys.argv[3]
OUT = sys.argv[4]
#SEQUENCES_PATH = 'gbff/Only_ORTHOLOGUES/' ## Path de la carpeta con los ortólogos.

output_file = '' ## Nombre del archivo de salida
#output_file = '.'.join(['Orthologues_Palindrome_sites.AllFrames',OUT,SPP,pattern,'txt'])
output_file = '.'.join(['Orthologues_Palindrome_sites.AllFrames',OUT,'txt'])
output = open (output_file, 'w') ## Abrimos el archivo de salida
output.write('FILE\tSpp\tSTART\tEND\tReadingFrame\tOrthLength\tPAL\tAA\n')

#output_fileM1 = '.'.join(['Orthologues_Palindrome_sites.M1',OUT,SPP,pattern,'txt'])
#outputM1 = open (output_fileM1, 'w') ## Abrimos el archivo de salida
#outputM1.write('FILE\tSpp\tSTART\tEND\tReadingFrame\tOrthLength\tPAL\tAA\n')

#output_fileM2 = '.'.join(['Orthologues_Palindrome_sites.M2',OUT,SPP,pattern,'txt'])
#outputM2 = open (output_fileM2, 'w') ## Abrimos el archivo de salida
#outputM2.write('FILE\tSpp\tSTART\tEND\tReadingFrame\tOrthLength\tPAL\tAA\n')

#output_fileM3 = '.'.join(['Orthologues_Palindrome_sites.M3',OUT,SPP,pattern,'txt'])
#outputM3 = open (output_fileM3, 'w') ## Abrimos el archivo de salida
#outputM3.write('FILE\tSpp\tSTART\tEND\tReadingFrame\tOrthLength\tPAL\tAA\n')

if SEQUENCES_PATH.endswith("/"): ## Esta parte la pongo por si corro este codigo en la terminal
    SequencesPath = re.sub('/', '', SEQUENCES_PATH) ## De este modo evito errores en el argumento path 
    SequencesDir = str(SEQUENCES_PATH)
else:
    SequencesPath = SEQUENCES_PATH
    SequencesDir = str("".join ([SEQUENCES_PATH,'/']))

def translate(seq):
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        '---':'-',
    }
    protein =""
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein+= table[codon]
    return protein

    
Orthologues = [x for x in os.listdir(SequencesDir) if x.endswith(".phy")] ## creo un arreglo con todos los ortólogis de la carpeta
#pattern = '[cC][-]*[gG][-]*[gG][-]*[cC][-]*[gG][-]*[cC][-]*[cC][-]*[gG]'
pattern = re.sub('G', '[gG][-]*', pattern)
pattern = re.sub('C', '[cC][-]*', pattern)
pattern = re.sub('T', '[tT][-]*', pattern)
pattern = re.sub('A', '[aA][-]*', pattern)
pattern = pattern[:-4]


k=0
RF1=0
RF2=0
RF3=0
for orthologue in Orthologues:
    file = re.sub('\.fna\.awk1\.mafft\.phy', '', orthologue)
    FNA = str("".join ([SequencesDir,orthologue]))
    spps = [str(record.description) for record in SeqIO.parse(open(FNA),'phylip')]
    sequencesDict = {str(record.description):str(record.seq) for record in SeqIO.parse(open(FNA),'phylip')}
    Sites = [match.span()[0] for match in re.finditer(pattern, sequencesDict[SPP])]
    EndSites = [match.span()[1] for match in re.finditer(pattern, sequencesDict[SPP])]
    j=0
    for site in Sites:
        k += 1
        for spp in spps:    
            end = EndSites[j]
            kmer = sequencesDict[spp][site:end]
            OrthLength = len(sequencesDict[spp])
            start = site
            mod1 = (start+1+2) % 3
            mod2 = (start+1+1) % 3
            mod3 = (start+1+0) % 3
            if mod1 == 0:
                RF = 1
                RF1+=1
            elif mod2 == 0:
                RF = 2
                RF2+=1
            elif mod3 == 0:
                RF = 3
                RF3+=1
                
            if len(kmer)==8 and RF ==1:
                AAstart = int(((start+3)/3)-1)
                AAend = AAstart+3
                AAseq = sequencesDict[spp]
                AA = translate(AAseq)[AAstart:AAend]
                codon1 = sequencesDict[spp][site:site+3]
                codon2 = sequencesDict[spp][site+3:site+3+3]
                codon3 = sequencesDict[spp][site+3+3:site+3+3+3]
                word = " ".join ([codon1,codon2,codon3])
                print ('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(file,spp,start+1,end,RF,OrthLength,word,AA))
                output.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(file,spp,start+1,end,RF,OrthLength,word,AA))
                #outputM1.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(file,spp,start+1,end,RF,OrthLength,word,AA))
            elif len(kmer)==8 and RF ==2:
                AAstart = int(((start+2)/3)-1)
                AAend = AAstart+3
                AAseq = sequencesDict[spp]
                AA = translate(AAseq)[AAstart:AAend]
                codon1 = sequencesDict[spp][site-1:site-1+3]
                codon2 = sequencesDict[spp][site-1+3:site-1+3+3]
                codon3 = sequencesDict[spp][site-1+3+3:site-1+3+3+3]
                word = " ".join ([codon1,codon2,codon3])
                print ('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(file,spp,start+1,end,RF,OrthLength,word,AA))
                output.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(file,spp,start+1,end,RF,OrthLength,word,AA))
                #outputM2.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(file,spp,start+1,end,RF,OrthLength,word,AA))
            elif len(kmer)==8 and RF ==3:
                AAstart = int(((start+1)/3)-1)
                AAend = AAstart+4
                AAseq = sequencesDict[spp]
                AA = translate(AAseq)[AAstart:AAend]
                codon1 = sequencesDict[spp][site-2:site-2+3]
                codon2 = sequencesDict[spp][site-2+3:site-2+3+3]
                codon3 = sequencesDict[spp][site-2+3+3:site-2+3+3+3]
                codon4 = sequencesDict[spp][site-2+3+3+3:site-2+3+3+3+3]
                word = " ".join ([codon1,codon2,codon3,codon4])
                print ('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(file,spp,start+1,end,RF,OrthLength,word,AA))
                output.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(file,spp,start+1,end,RF,OrthLength,word,AA))
                #outputM3.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(file,spp,start+1,end,RF,OrthLength,word,AA))
        j += 1
        print ('_________________________________________________________________________\n')
        #k = k + len(Sites)
output.close()
#outputM1.close()
#outputM2.close()
#outputM3.close()
print ("TOTAL: {} sitios".format(k))

print ("Hay {} sitios en el marco de lectura 1.".format(RF1/len(spps)))
print ("Hay {} sitios en el marco de lectura 2.".format(RF2/len(spps)))
print ("Hay {} sitios en el marco de lectura 3.".format(RF3/len(spps)))
