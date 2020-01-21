class Annotation(object):
    def __init__(self, gffLine, ranks):
        gffLine=gffLine.split("\t")
        self.chromosome=gffLine[ranks[0]]
        self.annotation=gffLine[ranks[2]]
        self.start=gffLine[ranks[3]]
        self.end=gffLine[ranks[4]]
        misc=gffLine[ranks[8]].strip("\n").split(";")
        self.misc={info.split("=")[0]:info.split("=")[1] for info in misc}
        self.id=self.misc["ID"]
        self.parents=self.misc["Parent"].split(",") if "Parent" in self.misc else None
        self.description=self.misc["description"] if "description" in self.misc else None

