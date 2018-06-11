from schrodinger.structutils.measure import measure_distance

def valid_sb(atom1, atom2):
    q1 = atom1.formal_charge
    q2 = atom2.formal_charge
    
    local_q1 = q1 + sum([n.formal_charge for n in atom1.bonded_atoms])
    local_q2 = q2 + sum([n.formal_charge for n in atom2.bonded_atoms])

    # (1) must have opposite formal charges
    # (2) formal charges must not be (approximately) cancelled by neighboring atoms
    return q1*q2 < 0 and local_q1*local_q2 < 0 and 'Zn' not in [atom1.element, atom2.element]

def valid_sb2(atom1, atom2):
    return atom1.formal_charge*atom2.formal_charge < 0

class SB:
    def __init__(self, atom1, atom2):
        # atom 1 and atom 2 must have opposite formal charges (see residue.py)
        self.res_atom = atom1
        self.lig_atom = atom2
        self.r = measure_distance(atom1, atom2) # atom1.dist_to(atom2)

    def score1(self):
        if self.r <= 3: return 1
        elif self.r <= 4.5: return (4.5 - self.r)/1.5
        return 0

    def score2(self):
        if self.r <= 4: return 1
        elif self.r <= 5: return (5 - self.r)#/1.5
        return 0

    def score3(self):
        if self.r <= 4.5: return 1
        elif self.r <= 6: return (6 - self.r)/1.5
        return 0

    def __str__(self):
        ra = self.res_atom
        la = self.lig_atom
        at1 = '\n+res_atom: {} {} \n++charge: {} \n++COO-: {}'.format(ra.index, ra.element, ra.formal_charge,0)#, is_carboxylate(ra)[0])
        at2 = '\n+lig_atom: {} {} \n++charge: {} \n++COO-: {}'.format(la.index, la.element, la.formal_charge,0)#, is_carboxylate(la)[0])
        return '\nSB score: {}, {} '.format(self.score(), self.newscore()) + at1 + at2 + '\n+r: ' + str(self.r)

class SB_Container:
    def __init__(self, lig, ind):
        self.lig = lig
        self.ind = ind
        self.all_sb = {}

    def add_residue(self, resnum, res):
        self.all_sb[resnum] = []
        for res_atom in res.chrg:
            for lig_atom in self.lig.chrg:
                if valid_sb2(res_atom, lig_atom):
                    self.all_sb[resnum].append(SB(res_atom, lig_atom))

    def filter_int(self):
        # enforces 1 sb per ligand formal charge
        used_fc = {}
        for r in self.all_sb:
            for sb in self.all_sb[r]:
                unique_index = (sb.lig_atom.index,r)
                if unique_index not in used_fc or used_fc[unique_index][0].r > sb.r:
                    used_fc[unique_index] = (sb, r)
        self.all_sb = {r:[sb for hid, (sb, sb_r) in used_fc.items() if sb_r == r] for r in self.all_sb}

    def score(self):
        all_scores = {}
        for r, sb_list in self.all_sb.items():
            for sb in sb_list:
                sc = [sb.score1(), sb.score2(), sb.score3()]
                for i,j in enumerate(self.ind):
                    all_scores[(j,r,'')] = all_scores.get(r, 0) + sc[i]
        return all_scores

    def __str__(self):
        return 'Salt Bridges: ' + ''.join([str(r) + '\n' + ''.join(['\n'+str(i) for i in self.all_sb[r]]) for r in sorted(self.all_sb.keys())])
