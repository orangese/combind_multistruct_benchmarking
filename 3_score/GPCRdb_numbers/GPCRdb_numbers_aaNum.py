import urllib,urllib2
import json
import re
from Bio.SeqUtils import seq1

""" function for mapping all GPCRdb numbers based on uniprot id """
class GPCRdbNumbers:
        def __init__(self, gpcrdb_generic_numbers_file = "/share/PI/rondror/docking/fingerprint/GPCRdb_numbers/All_species_gpcrdb_numbers_revised_17May2016.txt"):
                self.generic_numbers_dict = {}
                with open(gpcrdb_generic_numbers_file) as GENERIC:
                        for line in GENERIC:
                                (uniprot, aaNum, aaName, TM, generic_num) = line.rstrip().split("\t")
                                generic_num = re.sub("\.\d+", "", generic_num)
                                if uniprot in self.generic_numbers_dict:
                                        self.generic_numbers_dict[uniprot][aaNum] = (generic_num, aaName)
                                else:
                                        self.generic_numbers_dict[uniprot] = {}
                                        self.generic_numbers_dict[uniprot][aaNum] = (generic_num, aaName)
                self.aa_numbers_dict = self.create_aa_numbers()

        def create_aa_numbers(self):
                aa_numbers_dict = {}
                for uniprot in self.generic_numbers_dict:
                        aa_numbers_dict[uniprot] = {}
                        for aa_num in self.generic_numbers_dict[uniprot]:
                                generic_num, aaName = self.generic_numbers_dict[uniprot][aa_num]
                                aa_numbers_dict[uniprot][generic_num] = (aa_num, aaName)
                return aa_numbers_dict

        def contains(self, uniprot):
                return uniprot in self.generic_numbers_dict

        def get_generic_num(self, uniprot, aaNum):
            try:
                return self.generic_numbers_dict[uniprot][aaNum]
            except Exception:
                return None, None

        def get_aa_num(self, uniprot, generic_num):
                return self.aa_numbers_dict[uniprot][generic_num]
