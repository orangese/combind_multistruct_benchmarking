from score.density_estimate import DensityEstimate
from utils import mp
import pandas as pd

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

data_root = '/oak/stanford/groups/rondror/projects/ligand-docking/combind_vs/DUDE/combind'
mode = 'shape'

if mode == 'shape':
    stats_root = '/oak/stanford/groups/rondror/users/jpaggi/SHAPE_stats_new'
    template='{}/{}/scores/all_rd1_shape_{}/{}/shape_new.csv'
    out = '{}/{}_SHAPE_{}_{}.de'
    metric = 'SHAPE_mean'
else:
    assert mode == '2D'
    stats_root = '/oak/stanford/groups/rondror/users/jpaggi/2D_stats_new'
    template='{}/{}/scores/all_rd1_shape_{}/{}/similarity.csv'
    out = '{}/{}_2D_{}_{}.de'
    metric = 'active_sim_mean'


def compute_stats(protein, n):
    print('compute', protein)
    nat, ref = [], []
    for i in range(5):
        fname = template.format(data_root, protein, n, i)
        df = pd.read_csv(fname).set_index('ID')
        X = df[metric]
        active = df.index.str.contains('CHEMBL')
        _nat = DensityEstimate(domain=(0, 1), sd=0.01, reflect=True).fit(X[active])
        _ref = DensityEstimate(domain=(0, 1), sd=0.01, reflect=True).fit(X[~active])
        nat += [_nat]
        ref += [_ref]
    return protein, DensityEstimate.merge(nat), DensityEstimate.merge(ref)

for n in [5]: #[0, 1, 3, 5, 10]:
    stats = mp(compute_stats, [(x, n) for x in dude], 24)

    for protein in dude:
        print('merge', protein)
        nat = []
        ref = []
        for (_protein, _nat, _ref) in stats:
            if protein != _protein:
                nat += [_nat]
                ref += [_ref]
        nat = DensityEstimate.merge(nat)
        ref = DensityEstimate.merge(ref)

        nat.write(out.format(stats_root, 'nat', protein, n))
        ref.write(out.format(stats_root, 'ref', protein, n))