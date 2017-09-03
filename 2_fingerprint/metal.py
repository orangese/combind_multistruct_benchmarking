from schrodinger.structutils.measure import measure_distance

class Metal:
    def __init__(self, metal_atom, lig_atom):
        self.r = measure_distance(metal_atom, lig_atom)
        self.lig_atom = lig_atom

    def score(self):
        if self.r <= 2: return 1
        elif self.r <= 3: return (3 - self.r)
        else: return 0

    def __str__(self):
        str1 = '\nMetal score: ' + str(self.score())
        str3 = '\n dist: ' + str(self.r)
        str2 = '\nLigand atom: ' + self.lig_atom.element + ' ' + str(self.lig_atom.index)
        return str1+str3+str2+'\n'

class Metal_Container:
    def __init__(self, lig_st):
        self.lig_st = lig_st
        self.all_metal = {}

    def has_lone_pair(self, atom):
        if atom.element in ['O', 'S'] and atom.bond_total <= 2:
            return True
        elif atom.element == 'N' and atom.bond_total <= 3:
            return True
        return False

    def add_residue(self, resnum, res_st):
        for res_atom in res_st.atom:
            if not res_atom.element == 'Zn': continue
            if resnum not in self.all_metal: self.all_metal[resnum] = []
            for lig_atom in self.lig_st.atom:
                if not self.has_lone_pair(lig_atom): continue
                self.all_metal[resnum].append(Metal(res_atom, lig_atom))

    def filter_int(self):
        for r in self.all_metal:
            best_metal = None
            for metal in self.all_metal[r]:
                if best_metal == None or metal.score() > best_metal.score():
                    best_metal = metal
            self.all_metal[r] = [best_metal]
    
    def score(self):
        return {
            r : [
                sum([metal.score() for metal in self.all_metal[r]])
            ] for r in self.all_metal
        }
    def __str__(self):
        return 'Metal:\n' + ''.join([str(r) + ''.join(['\n'+str(i) for i in self.all_metal[r]]) for r in sorted(self.all_metal.keys())])

