class Annotation(object):
    """A parsed line of a GFF3 file."""
    def __init__(self, gffLine, ranks):
        """
        Arguments:
        gffLine (str) : Raw line of the GFF3 file
        ranks (list) : Order of fields as in GFF3 specs

        chromosome (str) : Chromosome name
        annotation (str) : Type of feature; 'gene', 'mRNA', 'CDS' and 'exon'
        start (int) : Start position (included) of feature
        end (int) : End position (included) of feature 
        strand (str) : Strand implicated among '+' and '-'
        phase (str) : Start of codon among '0', '1', '2' (CDS only) and '.'
        misc (dict) : GFF3 INFO splitted by category; 'ID=' is mandatory
        id (str) : Unique ID of the feature
        parents (str) : Parent(s) of the feature
        description (str) : Annotated function of the feature

        """
        gffLine=gffLine.split("\t")
        self.chromosome=gffLine[ranks[0]]
        self.annotation=gffLine[ranks[2]]
        self.start=int(gffLine[ranks[3]])
        self.end=int(gffLine[ranks[4]])
        self.strand=gffLine[ranks[6]]
        self.phase=gffLine[ranks[7]]
        misc=gffLine[ranks[8]].strip("\n").split(";")
        self.misc={info.split("=")[0]:info.split("=")[1] for info in misc}
        self.id=self.misc["ID"]
        self.parents=self.misc["Parent"].split(",") if "Parent" in self.misc else None
        self.description=self.misc["description"] if "description" in self.misc else None