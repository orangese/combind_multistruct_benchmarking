import urllib,urllib2
import json
import re
from Bio.SeqUtils import seq1

""" function for mapping all GPCRdb numbers based on uniprot id """
def Get_GPCRdb_Numbers():
	generic_numbers_dict = {}
	gpcrdb_generic_numbers_file = "./All_species_gpcrdb_numbers_revised_17May2016.txt"
	with open(gpcrdb_generic_numbers_file) as GENERIC:
		for line in GENERIC:
			(uniprot, aaNum, aaName, TM, generic_num) = line.rstrip().split("\t")
			generic_num = re.sub("\.\d+", "", generic_num)
			if uniprot in generic_numbers_dict.keys():
				generic_numbers_dict[uniprot][generic_num] = aaNum
			else:
				generic_numbers_dict[uniprot] = {}
				generic_numbers_dict[uniprot][generic_num] = aaNum
	return generic_numbers_dict
