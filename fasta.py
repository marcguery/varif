class Fasta(object):
    def __init__(self):
        self.data={}
        self.dnatable={ 
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
        }
        self.dnareverse={
            'A':'T', 'T':'A', 'C':'G', 'G':'C'
        }
    
    def load_data_from_FASTA(self, fasta):
        with open(fasta, 'r') as fastaFile:
            fastaInfo=[line.rstrip("\n") for line in fastaFile]
            fastaHeaders=[index for index, line in enumerate(fastaInfo) if line[0]==">"]
            n=1
            while n < len(fastaHeaders):
                currentIdx=fastaHeaders[n-1]
                nextIdx=fastaHeaders[n]
                chromosome=fastaInfo[currentIdx].split("|")[0].strip("> ")
                self.data[chromosome]="".join(fastaInfo[currentIdx+1:nextIdx])
                n+=1
            chromosome=fastaInfo[nextIdx].split("|")[0].strip("> ")
            self.data[chromosome]="".join(fastaInfo[nextIdx+1:])
    
    def reverse_strand(self, cds):
        cds=cds[::-1]
        newcds=""
        for i in range(len(cds)):
            newcds+=self.dnareverse[cds[i]]
        return newcds

    def translate_CDS(self, cds, strand, phase):
        protein=""
        aminoacid=""
        cds=self.reverse_strand(cds) if strand=="-" else cds
        i=phase
        n=len(cds)
        while i <= n-3 and aminoacid != "_":
            aminoacid=self.dnatable[cds[i:i+3]]
            protein+=aminoacid
            i+=3
        return protein



        