from .config import Config

class Families(object):
    """Store data loaded from a PED file enabling the possibility of grouping by family or lineage"""

    def __init__(self):
        """
        pedfile (str) : File path for the PED file
        families (dict) : Family name (key) and samples belonging to it (value)
        offspring (list) : Offspring IDs (same order as 'mates')
        mates (list) : Mates IDs (same order as 'offspring')
        samples (list) : All samples from the PED file
        
        """
        self.pedfile = None
        self.families = {}
        self.offspring = []
        self.mates = []
        self.samples = []
        
    
    def check_mates(self, clean = True):       
        """Checks the integrity of the mates and their corresponding offspring

        Args:
            clean (bool): Updates lists if an error is found
        """ 
        dellist = []
        for index, pair in enumerate(self.mates):
            for parent in pair:
                if parent not in self.samples:
                    Config.error_print("Parent %s not found in sample list"%(parent))
                    if index not in dellist:
                        dellist.append(index)
        
        if clean is True:
            for incr, index in enumerate(dellist):
                del self.mates[index-incr]
                del self.offspring[index-incr]
        
        if len(self.mates) != len(self.offspring):
            Config.error_print("Mates and offspring numbers are not matching")
            raise ValueError
            
    
    def build_families(self, sample, family):
        """Appends a sample to its corresponding family

        Args:
            sample (str): ID of the sample
            family (str): ID of the family
        """        
        if family in self.families:
            if sample not in self.families[family]:
                self.families[family].append(sample)
            else:
                Config.error_print("Sample %s from family %s is duplicated! Remove the duplicated lines or change the family/sample identifier"%(sample, family))
                raise ValueError()
        else:
            self.families[family] = [sample]
    
    def build_lineages(self, sample, father, mother):
        """Appends a sample to its corresponding lineage

        Args:
            sample (str): ID of the sample
            father (str): ID of the father
            mother (str): ID of the mother
        """        
        curr_mates = [min([father, mother]), max([father, mother])]
        
        mateid = None
        for index, pair in enumerate(self.mates):
            if curr_mates == pair:
                mateid = index
                break
        
        if mateid is not None:
            if sample not in self.offspring[mateid]:
                self.offspring[mateid].append(sample)
        else:
            self.mates.append(curr_mates)
            self.offspring.append([sample])
        
    
    def read_ped(self, ped):
        """
        Reads a PED file

        ped (str) : File path of PED file

        """
        with open(ped, 'r') as f:
            self.pedfile = f.readlines()
        for line in self.pedfile:
            data = line.split()
            family = data[0]
            sample = data[1]
            father = data[2]
            mother = data[3]
            #sex = data[4]
            #genotype = data[5]
            if family != "":
                self.build_families(sample, family)
            if father != "" or mother != "":
                self.build_lineages(sample, father, mother)
            self.samples.append(sample)
        
        self.check_mates()
