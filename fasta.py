class Fasta(object):
    def __init__(self):
        self.data={}
    
    def load_data_from_FASTA(self, fasta):
        with open(fasta, 'r') as fastaFile:
            fastaInfo=[line.rstrip("\n") for line in fastaFile]
            fastaHeaders=list(filter(lambda x:x[1][0]==">", enumerate(fastaInfo)))
            n=0
            while n < len(fastaHeaders):
                currentIdx=fastaHeaders[n-1][0]
                nextIdx=fastaHeaders[n][0]
                chromosome=fastaInfo[currentIdx].split("|")[0].strip("> ")
                self.data[chromosome]="".join(fastaInfo[currentIdx+1:nextIdx])
                n+=1
            chromosome=fastaInfo[nextIdx].split("|")[0].strip("> ")
            self.data[chromosome]="".join(fastaInfo[nextIdx+1:])


        