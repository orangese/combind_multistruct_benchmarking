from schrodinger.structutils.measure import measure_distance

def valid_sb(atom1, atom2):
    q1 = atom1.formal_charge
    q2 = atom2.formal_charge
    
    local_q1 = q1 + sum([n.formal_charge for n in atom1.bonded_atoms])
    local_q2 = q2 + sum([n.formal_charge for n in atom2.bonded_atoms])

    # (1) must have opposite formal charges
    # (2) formal charges must not be (approximately) cancelled by neighboring atoms
    return q1*q2 < 0 and local_q1*local_q2 < 0 and 'Zn' not in [atom1.element, atom2.element]

class SB:
    def __init__(self, atom1, atom2):
        # atom 1 and atom 2 must have opposite formal charges (see residue.py)
        self.res_atom = atom1
        self.lig_atom = atom2
        self.r = measure_distance(atom1, atom2) # atom1.dist_to(atom2)

    #def get_distance(self): # special treatment for carboxylate
    #    r1 = self.lig_atom.dist_to(self.res_atom)
    #    if is_carboxylate(self.res_atom)[0]:
    #        r2 = self.lig_atom.dist_to(is_carboxylate(self.res_atom)[1])
    #        return 2*r1*r2/(r1+r2) # effective distance
    #    elif is_carboxylate(self.lig_atom)[0]:
    #        r2 = self.res_atom.dist_to(is_carboxylate(self.lig_atom)[1])
    #        return 2*r1*r2/(r1+r2)
    #    return r1

    def score(self):
        # scales with 1/r (as in electric potential energy) and is capped at 2
        #if self.r <= 3: return 2
        #return 5 - self.r
        return min(2, 4.0/self.r)

    def newscore(self):
        if self.r <= 3: return 1
        elif self.r <= 4.5: return (4.5 - self.r)/1.5
        return 0

    def __str__(self):
        ra = self.res_atom
        la = self.lig_atom
        at1 = '\n+res_atom: {} {} \n++charge: {} \n++COO-: {}'.format(ra.index, ra.element, ra.formal_charge,0)#, is_carboxylate(ra)[0])
        at2 = '\n+lig_atom: {} {} \n++charge: {} \n++COO-: {}'.format(la.index, la.element, la.formal_charge,0)#, is_carboxylate(la)[0])
        return '\nSB score: {}, {} '.format(self.score(), self.newscore()) + at1 + at2 + '\n+r: ' + str(self.r)

class SB_Container:
    def __init__(self, lig_st, sub_st_map, ind):
        self.lig_st = lig_st
        self.sub_st_map = sub_st_map
        self.ind = ind
        self.all_sb = {}

    def add_residue(self, resnum, res_st):
        self.all_sb[resnum] = []
        for res_atom in res_st.atom:
            for lig_atom in self.lig_st.atom:
                if valid_sb(res_atom, lig_atom):
                    self.all_sb[resnum].append(SB(res_atom, lig_atom))

    def filter_int(self):
        # enforces 1 sb per ligand formal charge
        #pass        
        used_fc = {}
        for r in self.all_sb:
            for sb in self.all_sb[r]:
                unique_index = sb.lig_atom.index
                if unique_index not in used_fc or used_fc[unique_index][0].newscore() < sb.newscore():
                    used_fc[unique_index] = (sb, r)
        self.all_sb = {r:[sb for hid, (sb, sb_r) in used_fc.items() if sb_r == r] for r in self.all_sb}

    def score(self):
        all_scores = {}
        for r, sb_list in self.all_sb.items():
            for sb in sb_list:
                lig_sub_st = ','.join([str(j) for j in self.sub_st_map.get(sb.lig_atom.index, '')])
                key = (self.ind[0], r, lig_sub_st)
                all_scores[key] = all_scores.get(key, 0) + sb.newscore()
        return all_scores
#        return {
#            r : [
#                sum([sb.score() for sb in self.all_sb[r]]),
#                sum([sb.newscore() for sb in self.all_sb[r]])
#            ] for r in self.all_sb
#        }

    def __str__(self):
        return 'Salt Bridges: ' + ''.join([str(r) + '\n' + ''.join(['\n'+str(i) for i in self.all_sb[r]]) for r in sorted(self.all_sb.keys())])
