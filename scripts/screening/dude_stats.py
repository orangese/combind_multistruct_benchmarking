from score.statistics import *
from utils import mp, mkdir

dude = ["aa2ar","abl1","ace","aces","ada","ada17","adrb1","adrb2","akt1",
        "akt2","aldr","ampc","andr","aofb","bace1","braf","cah2","casp3","cdk2",
        "comt","cp2c9","cp3a4","csf1r","cxcr4","def","dhi1","dpp4","drd3","dyr",
        "egfr","esr1","esr2","fa10","fa7","fabp4","fak1","fgfr1","fkb1a","fnta",
        "fpps","gcr","glcm","gria2","grik1","hdac2","hdac8","hivint","hivpr",
        "hivrt","hmdh","hs90a","hxk4","igf1r","inha","ital","jak2","kif11",
        "kit","kith","kpcb","lck","lkha4","mapk2","mcr","met","mk01","mk10",
        "mk14","mmp13","mp2k1","nos1","nram","pa2ga","parp1","pde5a","pgh1",
        "pgh2","plk1","pnph","ppara","ppard","pparg","prgr","ptn1","pur2",
        "pygm","pyrd","reni","rock1","rxra","sahh","src","tgfr1","thb","thrb",
        "try1","tryb1","tysy","urok","vgfr2","wee1","xiap"]

bpp = ["5HT2B","AR","B1AR","B2AR","BACE1","BRD4","CDK2","DAT","DHFR","ELANE",
       "ERA","F10","F2","GLUT1","HSP90AA1","MGLUR5","NR3C1","NR3C2","P00760",
       "P19491","P22756","PDE10A","PLAU","PTPN1","PYGM","Q05586-Q12879",
       "SIGMAR1","SLC6A4","SMO","VDR"]

dude_to_bpp = {
    'andr': 'AR',
    'adrb1': 'B1AR',
    'adrb2': 'B2AR',
    'bace1': 'BACE1',
    'cdk2': 'CDK2',
    'dyr': 'DHFR',
    'esr1': 'ERA',
    'fa10': 'F10',
    'hs90a': 'HSP90AA1',
    'gcr': 'NR3C1',
    'mcr': 'NR3C2',
    'try1': 'P00760',
    'urok': 'PLAU',
    'ptn1': 'PTPN1',
    'grik1': 'P22756',
    'gria2': 'P19491',
    'pygm': 'PYGM',
}

for k, v in dude_to_bpp.items():
    assert k in dude, k
    assert v in bpp, v

data_root = '/oak/stanford/groups/rondror/users/jpaggi/combind2'
pairs_root = '/oak/stanford/groups/rondror/users/jpaggi/combind2_pairs'
stats_root = '/oak/stanford/groups/rondror/users/jpaggi/combind2_stats'
dude_stats_root = '/oak/stanford/groups/rondror/users/jpaggi/dude_stats'
features = ['shape', 'mcss', 'hbond', 'saltbridge', 'contact']

args = []
for protein in bpp:
    args += [(protein, data_root, pairs_root)]
mp(pair_features, args, 12)

args = []
for protein in bpp:
    mkdir('{}/{}'.format(stats_root, protein))
    args += [(protein, pairs_root, stats_root, features)]
mp(compute_stats, args, 12)

args = []
for protein in dude:
    proteins = [_protein for _protein in bpp
                if _protein not in dude_to_bpp or _protein != dude_to_bpp[protein]]

    root = '{}/{}'.format(dude_stats_root, protein)
    mkdir(root)
    fname = root + '/{}_{}.de'
    args += [(proteins, stats_root, fname, features)]
mp(merge_stats, args, 12)
