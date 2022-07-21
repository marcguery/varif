from sys import stderr

class Fasta(object):
    """Handle FASTA files and basic operations on bases."""
    def __init__(self):
        """
        data (dict) : Key (chromosome), value (sequence)
        dnatable (dict) : Map of codons and aminoacids
        dnareverse (dict) : Map of complementary bases

        """
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
        """
        Store a list of chromosomes from FASTA

        fasta (str) : File path of FASTA file

        """
        with open(fasta, 'r') as fastaFile:
            fastaInfo=[line.rstrip("\n") for line in fastaFile]
            #Looking for headers
            fastaHeaders=[index for index, line in enumerate(fastaInfo) if line[0]==">"]
            #Iterating over headers only
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
        """
        Get the complement of a DNA sequence

        cds (str) : Upper case of base sequences

        return (str) : The reverted sequence
        """
        cds=cds[::-1]
        newcds=""
        for i in range(len(cds)):
            newcds+=self.dnareverse[cds[i]] if cds[i] in self.dnareverse else '?'
        return newcds
    
    def insert_mutation(self, chromosome, genesequence, strand, startgene, startref, reference, mutation):
        """
        Update a CDS sequence by replacing a reference sequence by a mutation

        chromosome (str) : Name of the chromosome where the CDS is
        genesequence (str) : Upper case sequence of the CDS
        strand (str) : Either "+" or "-" indicating the direction of translation
        startgene (int) : Chromosomal position of the first base of the CDS
        startref (int) : Chromosomal position of the first base of the reference sequence
        reference (str) : Upper case sequence of the reference sequence
        mutation (str) : Upper case sequence of the mutation

        return (list) : Sequence of the new CDS, relative positions of the mutation

        """
        end=startgene + len(genesequence)
        posIndex=startref - startgene
        startIndex = max(0, posIndex)
        endIndex=min(len(genesequence), startIndex+len(reference))
        #Add bases if the mutation is shorter than the reference
        shift = len(reference) - len(mutation)
        #Hanging reference bases from the reference are not included in the size difference
        #Ref bases before the start of CDS
        upshift = shift + posIndex +1 if posIndex < -1 else shift
        #Ref bases after the end of CDS
        downshift = shift - (posIndex+len(reference)-len(genesequence)) if posIndex + len(reference) >= len(genesequence) else shift

        #Remove bases from the mutation if they increase the initial CDS size by adding bases before the first codon
        # or add bases from the up/downstream sequence to match the initial CDS size
        #'+' strand: remove mutated bases if located before the start of the CDS
        upmutshift = max(0, -upshift) if startref <= startgene and strand == "+" else 0
        #'+' strand: add upstream bases if mutation is too short and reference starts before the start of CDS
        # '-' strand: add upstream bases if mutation is shorter than the reference
        upshift = max(0, upshift) if strand == "-" and startIndex + len(mutation) < len(genesequence) or startref <= startgene else 0
        #'-' strand: remove mutated bases if located after the end of the CDS
        downmutshift = max(0, -downshift) if startIndex + len(mutation) >= len(genesequence) and strand == "-" else 0
        #'-' strand: add upstream bases if mutation is too short and reference ends after the end of the CDS
        # '+' strand: add downstream bases if mutation is shorter than the reference
        downshift = max(0, downshift) if strand == "+" and startref > startgene or startIndex + len(reference) >= len(genesequence) else 0
            
        genesequencemut = self.data[chromosome][startgene-upshift:startgene]
        genesequencemut+= genesequence[0:startIndex]
        genesequencemut+= mutation[upmutshift:len(mutation)-downmutshift]
        genesequencemut+= genesequence[endIndex:]
        genesequencemut+= self.data[chromosome][end:end+downshift]

        #Add more bases to make the CDS a multiple of 3
        # from upstream (- strand) or downstream (+ strand)
        basestoaddforcds = -len(genesequencemut)%3
        genesequencemut = self.data[chromosome][startgene-upshift-basestoaddforcds:startgene-upshift]+genesequencemut if strand == "-" else genesequencemut+self.data[chromosome][end+downshift:end+downshift+basestoaddforcds]

        startIndexMut = startIndex+upshift+basestoaddforcds if strand == "-" else startIndex+upshift
        endIndexMut=min(len(genesequencemut), startIndexMut+len(mutation[upmutshift:len(mutation)-downmutshift]))
        return [genesequencemut, startIndexMut, endIndexMut]
    
    def window_sequence(self, chromosome, position, sequence, windowBefore, windowAfter):
        """
        Extract DNA bases before and after a given sequence

        chromosome (str) : Name of the chromosome
        position (int) : Chromosomal position of the first base of the sequence
        sequence (str) : Upper case sequence of the sequence
        windowBefore (int) : Number of bases to extract before the sequence
        windowAfter (int) : Number of bases to extract after the sequence

        return (list) : Bases before and after the sequence and the sequence itself

        """
        seqLength = len(sequence)
        posBefore = position - windowBefore
        posAfter = position + windowAfter
        if posBefore < 1 or posAfter + seqLength > len(self.data[chromosome]):
            print("DNA window too wide at position (%s:%s). Changing to maximal window..."%(chromosome, position), file=stderr)
            posBefore = 1 if posBefore < 1 else posBefore
            posAfter = -2 - seqLength if posAfter + seqLength > len(self.data[chromosome]) else posAfter
        return [self.data[chromosome][posBefore-1:position-1],self.data[chromosome][position+seqLength-1:posAfter+seqLength-1]]
    
    def window_CDS(self, genesequence, strand, startIndex, endIndex, windowBefore, windowAfter):
        """
        Extract codons before and after initial codon(s) sharing at least one base of a given sequence

        chromosome (str) : Name of the chromosome
        genesequence (str) : Upper case sequence of the CDS
        strand (str) : Either "+" or "-" indicating the direction of translation
        startIndex (int) : Relative position of the first base of the sequence
        endIndex (int) : Relative position of the last base of the sequence
        windowBefore (int) : Number of codons to extract before the intiial codon
        windowAfter (int) : Number of codons to extract after the intiial codon

        return (list) : Codons, the CDS position of the first initial codon and the phase

        """
        #The number of nucleotides to skip
        # knowing that the window is +/- windows nucleotides after mutation
        #Translation starts at start of feature
        if strand=="+":
            shiftBefore=min(startIndex, windowBefore)
            shiftAfter=min(len(genesequence)-endIndex, windowAfter)
            aaPos=int(startIndex/3)+1
            phase=-(startIndex-shiftBefore)%3
        #Translation starts at end of feature
        elif strand=="-":
            shiftBefore=min(startIndex, windowAfter)
            shiftAfter=min(len(genesequence)-endIndex, windowBefore)
            aaPos=int((len(genesequence)-endIndex)/3)+1
            phase=-(len(genesequence)-(endIndex+shiftAfter))%3
        
        cds=genesequence[startIndex-shiftBefore:endIndex+shiftAfter]
        return [cds, aaPos, phase]

    def translate_CDS(self, cds, strand, phase):
        """
        Translate a DNA strand to a protein

        cds (str) : The upper case list of bases
        strand (str) : The strand to tell how to read CDS; '+' or '-'
        phase (int) : The number of bases to skip before the first codon

        return (str) : The upper case list of aminoacids
        """
        protein=""
        aminoacid=""
        cds=self.reverse_strand(cds) if strand=="-" else cds
        i=phase
        n=len(cds)
        while i <= n-3 and aminoacid != "_":
            aminoacid=self.dnatable[cds[i:i+3]] if cds[i:i+3] in self.dnatable else '?'
            protein+=aminoacid
            i+=3
        return protein
