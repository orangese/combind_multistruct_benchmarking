import os
from schrodinger.structure import StructureReader
from schrodinger.structutils.rmsd import renumber_conformer, calculate_in_place_rmsd
import centroid

def get_alignment_stats():
    check_ligand_alignment()
    check_binding_pocket_alignment()

def check_ligand_alignment():
    centroids = {}
    for f in os.listdir('aligned_ligands'):
        lig = StructureReader('aligned_ligands/{}'.format(f)).next()
        centroids[f.split('.')[0]] = centroid.average_atom_pos(lig)

    total_num_ligs = len(centroids.keys())

    def get_worst(centroids):
        bad_pairs = {}
        for i in range(len(centroids.keys())):
            for j in range(i+1, len(centroids.keys())):
                l1 = centroids.keys()[i]
                l2 = centroids.keys()[j]
                d = centroid.dist(centroids[l1], centroids[l2])
                if d > 8.0:
                    bad_pairs[l1] = bad_pairs.get(l1, 0) + 1
                    bad_pairs[l2] = bad_pairs.get(l2, 0) + 1
        if bad_pairs == {}: return None
        else:
            worst = sorted(bad_pairs.keys(), key=lambda x:-bad_pairs[x])[0]
            print worst, bad_pairs[worst]
            if bad_pairs[worst] == 1: return None
            return worst

    exclude = get_worst(centroids)
    with open('excluded_ligands.txt', 'w') as f:
        while exclude is not None:
            f.write(exclude+'\n')
            del centroids[exclude]
            exclude = get_worst(centroids)

    centroid_mean = [str(sum([centroids[l][i] for l in centroids.keys()])/len(centroids.keys())) for i in range(3)]
    with open('grid_center.txt', 'w') as f:
        f.write(','.join(centroid_mean))

def check_binding_pocket_alignment():
    residues = {}
    structs = {}
    for f in os.listdir('aligned_proteins'):
        name = f.split('.')[0]
        structs[name] = StructureReader('aligned_proteins/{}'.format(f)).next()
        for r in structs[name].residue:
            if r.hasMissingAtoms(): continue
            if (r.resnum, r.pdbres) not in residues:
                residues[(r.resnum, r.pdbres)] = {}
            residues[(r.resnum, r.pdbres)][name] = r.extractStructure()

    with open('residue_alignments.txt', 'w') as f:
        for r in residues:
            for i1 in range(len(structs.keys())):
                if structs.keys()[i1] not in residues[r]: continue
                for i2 in range(i1+1,len(structs.keys())):
                    if structs.keys()[i2] not in residues[r]: continue
                    r1 = residues[r][structs.keys()[i1]]
                    r2 = residues[r][structs.keys()[i2]]
                    if len(r1.atom) != len(r2.atom): continue
                    renumber_conformer(r1, r2, use_symmetry=True)
                    a1 = [a.index for a in r1.atom]
                    a2 = [a.index for a in r2.atom]
                    rmsd = calculate_in_place_rmsd(r1, a1, r2, a2, use_symmetry=True)
                    f.write('{} {} {} {} {}\n'.format(r[0],r[1], structs.keys()[i1], structs.keys()[i2], rmsd))
