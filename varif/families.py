class Families(object):
    """Store data loaded from a PED file enabling the possibility of grouping by family or lineage"""

    def __init__(self):
        """
        pedfile (str) : File path for the PED file
        families (dict) : Family name (key) and samples belonging to it (value)
        parents (dict) : Sample names for father and mother concatenated with blank space (key) and their offspring (value)
        sample (list) : All samples from the PED file
        
        """
        self.pedfile=None
        self.families={}
        self.parents={}
        self.samples=[]
    
    def read_ped(self, ped):
        """
        Reads a PED file

        ped (str) : File path of PED file

        """
        with open(ped, 'r') as f:
            self.pedfile=f.readlines()
        for line in self.pedfile:
            data=line.split()
            family=data[0]
            sample=data[1]
            father=data[2] if data[2] != "0" else ""
            mother=data[3] if data[3] != "0" else ""
            if father == "" and mother == "":
                parent = "NA"
            else:
                parent=min([father, mother])+" "+max([father, mother])
            sex=data[4]
            genotype=data[5]
            if family in self.families:
                if sample not in self.families[family]:
                    self.families[family].append(sample)
                else:
                    raise ValueError("Sample %s from family %s is duplicated! Remove the duplicated lines or change the family/sample identifier"%(sample, family))
            else:
                self.families[family] = [sample]
            if parent in self.parents:
                if sample not in self.parents[parent]:
                    self.parents[parent].append(sample)
            else:
                self.parents[parent] = [sample]
            self.samples.append(sample)
