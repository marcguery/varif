from .config import Config

class Fasta(object):
    """Handle FASTA files and basic operations on bases."""
    def __init__(self):
        """
        data (dict) : Key (chromosome), value (sequence)
        dnatable (dict) : Map of codons and aminoacids
        dnareverse (dict) : Map of complementary bases

        """
        self.data = {}
        self.dnatable = { 
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
        self.dnareverse = {
            'A':'T', 'T':'A', 'C':'G', 'G':'C'
        }
    
    def load_data_from_FASTA(self, fasta):
        """
        Store a list of chromosomes from FASTA

        fasta (str) : File path of FASTA file

        """
        with open(fasta, 'r') as fastaFile:
            fastaInfo = [line.rstrip("\n") for line in fastaFile]
            #Looking for headers
            fastaHeaders = [index for index, line in enumerate(fastaInfo) if line[0] == ">"]
            #Iterating over headers only
            n = 1
            while n < len(fastaHeaders):
                currentIdx = fastaHeaders[n-1]
                nextIdx = fastaHeaders[n]
                chromosome = fastaInfo[currentIdx].split("|")[0].strip("> ")
                self.data[chromosome] = "".join(fastaInfo[currentIdx+1:nextIdx])
                n += 1
            chromosome = fastaInfo[nextIdx].split("|")[0].strip("> ")
            self.data[chromosome] = "".join(fastaInfo[nextIdx+1:])
    
    def reverse_strand(self, cds):
        """
        Get the complement of a DNA sequence

        cds (str) : Upper case of base sequences

        return (str) : The reverted sequence
        
        """
        cds = cds[::-1]
        newcds = ""
        for i in range(len(cds)):
            newcds += self.dnareverse[cds[i]] if cds[i] in self.dnareverse else '?'
        return newcds
    
    def get_relative_gene_position(self, cdsCoords, position):
        """
        Transform a chromosomal position to a gene position (without intronic sequences)

        cdsCoords (list) : 0-based coordinate ranges of all the CDS sequences (ordered)
        position (int) : 0-based position to convert

        return (int) : Relative position (negative or positive) to the 0-based gene starting position 

        """
        indexCoord = 0
        sumCDS = 0
        if position < cdsCoords[indexCoord][0]:
            relposition = position - cdsCoords[indexCoord][0]
            return relposition
        while position not in range(cdsCoords[indexCoord][0], cdsCoords[indexCoord][1]):
            sumCDS += cdsCoords[indexCoord][1] - cdsCoords[indexCoord][0]
            indexCoord += 1
            if indexCoord == len(cdsCoords):
                if position < cdsCoords[len(cdsCoords)-1][0]:
                    return None
                else:
                    break
        relposition = sumCDS + position - cdsCoords[indexCoord][0]
        return relposition        
    
    def get_introns_coords(self, cdsCoords):
        """
        Retrieve the range coordinates of intronic sequences given CDS coordinates

        cdsCoords (list) : 0-based coordinate ranges of the CDS sequences

        return (list) : 0-based coordinate ranges of the intronic sequences

        """
        intronsCoords = []
        for indexcds in range(1, len(cdsCoords)):
            endpreviouscds = cdsCoords[indexcds-1][1]
            startcds = cdsCoords[indexcds][0]
            intronsCoords.append([endpreviouscds, startcds])
        return intronsCoords
            
    def remove_bases_from_features(self, sequence, featuresCoord, sequencePosition):
        """
        Remove bases from a sequence if they are located within given coordinate ranges

        sequence (str) : The sequence to remove bases from
        featuresCoord (list) : 0-based coordinate ranges of the features
        sequencePosition (int) : 0-based start position of the sequence

        return (list) : Sequence without bases located within features and new starting position of the sequence

        """
        coordtoremove = []
        sequenceLastposition = sequencePosition + len(sequence)
        newSequencePosition = sequencePosition
        indexSequence = 0
        indexFeature = 0

        #Will save from left to right the coordinates of the sequence that
        # are included in the range coordinates features
        while indexFeature < len(featuresCoord):
            currentPosition = max(indexSequence + sequencePosition, featuresCoord[indexFeature][0])
            if currentPosition >= featuresCoord[indexFeature][1]:
                #Sequence to check is after last base of feature
                indexFeature += 1
                continue
            if currentPosition >= len(sequence)+sequencePosition or sequenceLastposition <= featuresCoord[indexFeature][0]:
                #Sequence was completely checked
                break
            if currentPosition >= featuresCoord[indexFeature][0] and currentPosition < featuresCoord[indexFeature][1]:
                #Sequence to check is included in a feature
                #Store index of sequence located within the feature
                maxPosition = min(sequenceLastposition - sequencePosition, featuresCoord[indexFeature][1] - sequencePosition)
                coordtoremove.append([currentPosition - sequencePosition, maxPosition])
                if currentPosition == sequencePosition:
                    #If the firs bases of the sequence is removed, 
                    # shift the start position
                    newSequencePosition = sequencePosition + maxPosition
                indexFeature += 1
                indexSequence = maxPosition
        
        if coordtoremove == []:
            #No bases included in any feature
            return [sequence, sequencePosition]
        #Create a new sequence with the complentary intervals
        complementaryInterval = [0,coordtoremove[0][0]]
        newSequence = sequence[complementaryInterval[0]:complementaryInterval[1]]

        for indexcoord in range(1,len(coordtoremove)):
            complementaryInterval = [coordtoremove[indexcoord-1][1],coordtoremove[indexcoord][0]]
            newSequence += sequence[complementaryInterval[0]:complementaryInterval[1]]
        
        complementaryInterval = [coordtoremove[len(coordtoremove)-1][1],len(sequence)]
        newSequence += sequence[complementaryInterval[0]:complementaryInterval[1]]

        return [newSequence, newSequencePosition]
    
    def insert_mutation(self, chromosome, genesequence, strand, startgene, endgene, startref, reference, mutation):
        """
        Update a CDS sequence by replacing a reference sequence by a mutation

        chromosome (str) : Name of the chromosome where the merged CDS sequence is
        genesequence (str) : Upper case of the merged CDS sequence
        strand (str) : Either "+" or "-" indicating the direction of translation
        startgene (int) : Chromosomal position of the first base of the merged CDS sequence
        startref (int) : Relative position of the reference sequence to the 0-based starting position of the merged CDS sequence
        reference (str) : Upper case sequence of the reference sequence
        mutation (str) : Upper case sequence of the mutation

        return (list) : Sequence of the new merged CDS sequence, relative positions of the mutation

        """
        startIndex = max(0, startref)
        endIndex = min(len(genesequence), startIndex+len(reference))
        #Add bases if the mutation is shorter than the reference
        shift = len(reference) - len(mutation)
        #Hanging reference bases from the reference are not included in the size difference
        #Ref bases before the start of CDS
        upshift = shift + startref +1 if startref < -1 else shift
        #Ref bases after the end of CDS
        downshift = shift - (startref+len(reference)-len(genesequence)) if startref + len(reference) >= len(genesequence) else shift

        #Remove bases from the mutation if they increase the initial CDS size by adding bases before the first codon
        # or add bases from the up/downstream sequence to match the initial CDS size
        #'+' strand: remove mutated bases if located before the start of the CDS
        upmutshift = max(0, -upshift) if startref <= 0 and strand == "+" else 0
        #'+' strand: add upstream bases if mutation is too short and reference starts before the start of CDS
        # '-' strand: add upstream bases if mutation is shorter than the reference
        upshift = max(0, upshift) if strand == "-" and startIndex + len(mutation) < len(genesequence) or startref <= 0 else 0
        #'-' strand: remove mutated bases if located after the end of the CDS
        downmutshift = max(0, -downshift) if startIndex + len(mutation) >= len(genesequence) and strand == "-" else 0
        #'-' strand: add upstream bases if mutation is too short and reference ends after the end of the CDS
        # '+' strand: add downstream bases if mutation is shorter than the reference
        downshift = max(0, downshift) if strand == "+" and startref > 0 or startIndex + len(reference) >= len(genesequence) else 0
            
        genesequencemut = self.data[chromosome][startgene-upshift:startgene]
        genesequencemut+= genesequence[0:startIndex]
        genesequencemut+= mutation[upmutshift:len(mutation)-downmutshift]
        genesequencemut+= genesequence[endIndex:]
        genesequencemut+= self.data[chromosome][endgene:endgene+downshift]

        #Add more bases to make the CDS a multiple of 3
        # from upstream (- strand) or downstream (+ strand)
        basestoaddforcds = -len(genesequencemut)%3
        genesequencemut = self.data[chromosome][startgene-upshift-basestoaddforcds:startgene-upshift]+genesequencemut if strand == "-" else genesequencemut+self.data[chromosome][endgene+downshift:endgene+downshift+basestoaddforcds]

        startIndexMut = startIndex+upshift+basestoaddforcds if strand == "-" else startIndex+upshift
        endIndexMut = min(len(genesequencemut), startIndexMut+len(mutation[upmutshift:len(mutation)-downmutshift]))
        return [genesequencemut, startIndexMut, endIndexMut]
    
    def merge_CDS(self, chromosome, cdsCoords):
        """
        Create the spliced CDS sequence from a genomic CDS sequence

        chromosome (str) : Name of the chromosome
        cdsCoords (list) : 0-based coordinate ranges of CDS sequences (ordered)

        return (str) : Sequence of the merged CDS

        """
        
        genesequence = ""
        for coord in cdsCoords:
            genesequence += self.data[chromosome][coord[0]:coord[1]]
        assert len(genesequence)%3 == 0
        return genesequence
    
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
            Config.error_print("DNA window too wide at position (%s:%s). Changing to maximal window..."%(chromosome, position))
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
        windowBefore (int) : Number of codons to extract before the initial codon
        windowAfter (int) : Number of codons to extract after the initial codon

        return (list) : Codons, the CDS position of the first initial codon and the phase

        """
        #The number of nucleotides to skip
        # knowing that the window is +/- windows nucleotides after mutation
        #Translation starts at start of feature
        if strand == "+":
            shiftBefore = min(startIndex, windowBefore)
            shiftAfter = min(len(genesequence)-endIndex, windowAfter)
            aaPos = int(startIndex/3)+1
            phase = -(startIndex-shiftBefore)%3
        #Translation starts at end of feature
        elif strand == "-":
            shiftBefore = min(startIndex, windowAfter)
            shiftAfter = min(len(genesequence)-endIndex, windowBefore)
            aaPos = int((len(genesequence)-endIndex)/3)+1
            phase = -(len(genesequence)-(endIndex+shiftAfter))%3
        
        cds = genesequence[startIndex-shiftBefore:endIndex+shiftAfter]
        return [cds, aaPos, phase]

    def translate_CDS(self, cds, strand, phase):
        """
        Translate a DNA strand to a protein

        cds (str) : The upper case list of bases
        strand (str) : The strand to tell how to read CDS; '+' or '-'
        phase (int) : The number of bases to skip before the first codon

        return (str) : The upper case list of aminoacids
        
        """
        protein = ""
        aminoacid = ""
        cds = self.reverse_strand(cds) if strand == "-" else cds
        i = phase
        n = len(cds)
        while i <= n-3 and aminoacid != "_":
            aminoacid = self.dnatable[cds[i:i+3]] if cds[i:i+3] in self.dnatable else '?'
            protein += aminoacid
            i += 3
        return protein
