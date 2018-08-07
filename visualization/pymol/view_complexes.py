from pymol import cmd


def load_complexes(protein):
    cmd.delete('*')
    with open('bpp_outputs/structs.tsv') as fp:
        for line in fp:
            try:
                prot, grid, lig = line.strip().split('\t')
                if protein == prot:
                    cmd.load('bpp_data/{}/structures/proteins/{}_prot.mae'.format(prot, lig))
                    cmd.load('bpp_data/{}/structures/ligands/{}_lig.mae'.format(prot, lig))
                    if grid == lig:
                        cmd.util.cbao("{}_prot".format(lig))
                        cmd.util.cbay("het and {}_lig".format(lig))
                        cmd.show('sticks', "het and {}_lig".format(lig))
                    else:
                        cmd.util.cbab("{}_prot".format(lig))
                        cmd.util.cbag("het and {}_lig".format(lig))

                    cmd.show('spheres', 'het and {}_prot and ({}_lig expand 5)'.format(lig, lig))
            except:
                print line.strip()

    cmd.show('cartoon')
    cmd.set('cartoon_oval_length', '0.5')
    cmd.set('cartoon_transparency', '0.5')
    cmd.hide('everything', 'element H and not (element N+O extend 1)')
    cmd.hide('everything', 'name H')

cmd.extend('load_complexes', load_complexes)
