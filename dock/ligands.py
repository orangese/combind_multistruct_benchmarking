import os
import subprocess
from schrodinger.structure import StructureReader, StructureWriter

def process(input_file, output_file):
    with StructureReader(input_file) as reader, \
        StructureWriter(output_file) as writer:
        for st in reader:
            # Remove explicit stereochemistry specifications. These cause
            # errors in downstream steps.
            for k in st.property.keys():
                if 's_st_EZ_' in k or 'Chiral' in k:
                    st.property.pop(k)

            # Give each atom a unique name, ligands generated from smiles
            # strings will not have any atom name by default.
            names = set()
            counts = {}
            for atom in st.atom:
                if not atom.pdbname.strip():
                    if atom.element not in counts: counts[atom.element] = 0
                    counts[atom.element] += 1
                    atom.pdbname = atom.element + str(counts[atom.element])
                    
                    assert atom.pdbname not in names, atom.pdbname
                    names.add(atom.pdbname)
            writer.append(st)

def prep_ligand(root, name, smiles):
    smi_file = '{}/{}.smi'.format(root, name)
    mae_noname_file = '{}/{}_nonames.mae'.format(root, name)
    mae_file = '{}/{}.mae'.format(root, name)
    cmd = 'ligprep -WAIT -epik -ismi {} -omae {}'.format(smi_file, mae_noname_file)

    with open(smi_file, 'w') as fp:
        fp.write('{} {}\n'.format(smiles, name))

    subprocess.run(cmd, shell=True, cwd=root)
    process(mae_noname_file, mae_file)

def prep_ligands(lm):
    unfinished = []
    for name, info in lm.pdb.items():
        if name not in lm.prepped:
            unfinished += [(name, info['SMILES'])]

    root = lm.path('LIGANDS_ROOT')

    if unfinished:
        print('Processing {} ligands'.format(len(unfinished)))
        os.system('mkdir -p {}'.format(root))
        
        for name, smiles in unfinished:
            path = '{}/{}'.format(root, name)
            os.system('rm -rf {}'.format(path))
            os.system('mkdir {}'.format(path))

            prep_ligand(path, name, smiles)
