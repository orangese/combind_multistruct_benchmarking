from pymol import cmd
import sys
import os
from glob import glob

def load_results(directory):
    cmd.delete('*')
    for f in glob("{}/mcss_*".format(directory)):
        cmd.load(f)
    
    cmd.load("{}/ligs.mae".format(directory))
    cmd.hide('everything', 'element H')
    cmd.util.cbab("ligs")

    with open("{}/mcss.txt".format(directory)) as fp:
        for line in fp:
            print line.strip()

cmd.extend('load_results', load_results)
