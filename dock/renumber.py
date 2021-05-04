from schrodinger.structutils.transform import get_centroid as _get_centroid
import numpy as np

def hetalign(ref, query):
    """
    Assign het groups with identical chemical structure that are positioned
    similarly the same chain, residue name, residue number, and atom names.
    """
    return query

################################################################################

def score_matrix(S, gap_cost=1):
    H = np.zeros((S.shape[0]+1, S.shape[1]+1))
    T = -1+np.zeros(H.shape)

    for i in range(1, H.shape[0]):
        for j in range(1, H.shape[1]):
            match = H[i-1, j-1] + S[i-1, j-1]
            delete = H[i-1, j] - gap_cost
            insert = H[i, j-1] - gap_cost

            if match > max(delete, insert, 0):
                H[i, j] = match
                T[i, j] = 0
            elif delete > max(insert, 0):
                H[i, j] = delete
                T[i, j] = 1
            elif insert > 0:
                H[i, j] = insert
                T[i, j] = 2
    return H, T

def traceback(H, T, initial=True):
    """
    Return a pair of lists the length of the alignment. In the ith position,
    place the index of the sequence aligned at that position or None if there
    is a gap.
    """
    if H[-1, -1] == 0:
        return [], []

    if initial:
        # Find the argmax of H. This is where we will start the traceback
        i, j = np.unravel_index(H.argmax(), H.shape)
    else:
        if T[-1, -1] == 0:
            i, j = T.shape[0]-1, T.shape[1]-1
        elif T[-1, -1] == 1:
            i, j = T.shape[0]-1, None
        elif T[-1, -1] == 2:
            i, j = None, T.shape[1]-1
        else:
            assert False
    remaining = traceback(H[:i, :j], T[:i, :j], initial=False)
    return remaining[0]+[i], remaining[1]+[j]

################################################################################

def distance(res1, res2):
    dist = np.sqrt(np.sum((res1[-1] - res2[-1])**2))
    score = np.exp(5-dist) / (1 + np.exp(5-dist))
    return 8*score - 4

def distance_matrix(chain1, chain2):
    distances = np.zeros((len(chain1), len(chain2)))
    for i, res1 in enumerate(chain1):
        for j, res2 in enumerate(chain2):
            distances[i, j] = distance(res1, res2)
    return distances

def substitution(res1, res2):
    match = int(res1[-2] == res2[-2])
    return 6*match - 3

def substitution_matrix(chain1, chain2):
    substitutions = np.zeros((len(chain1), len(chain2)))
    for i, res1 in enumerate(chain1):
        for j, res2 in enumerate(chain2):
            substitutions[i, j] = substitution(res1, res2)
    return substitutions

################################################################################

def standardize_resname(s):
    s = s.strip()
    resnames = {
        'ALA': 'A',
        'ARG': 'R',
        'ASN': 'N',
        'CYS': 'C',
        'GLN': 'Q',
        'GLY': 'G',
        'ILE': 'I',
        'LEU': 'L',
        'MET': 'M',
        'PHE': 'F',
        'PRO': 'P',
        'SER': 'S',
        'THR': 'T',
        'TRP': 'W',
        'TYR': 'Y',
        'VAL': 'V',

        'ASP': 'D',
        'ASH': 'D',
        'GLU': 'E',
        'GLH': 'E',
        'HIS': 'H',
        'HSD': 'H',
        'HID': 'H',
        'HSE': 'H',
        'HIE': 'H',
        'HSP': 'H',
        'HIP': 'H',
        'LYS': 'K',
        'LYN': 'K',
    }

    if s in resnames:
        return resnames[s]
    return s

def get_centroid(residue):
    st = residue.extractStructure()
    return _get_centroid(st)[:3]

def get_residues(st, centroid=True):
    residues = []
    for residue in st.residue:
        if centroid:
            residues += [(residue.chain,
                          residue.resnum,
                          residue.inscode,
                          standardize_resname(residue.pdbres),
                          get_centroid(residue))
                        ]
        else:
            residues += [(residue.chain,
                          residue.resnum,
                          residue.inscode,
                          standardize_resname(residue.pdbres))
                        ]
    return sorted(residues)

################################################################################

def align_chains(ref_chain, query_chain):
    ref_res = get_residues(ref_chain)
    query_res = get_residues(query_chain)

    d = distance_matrix(ref_res, query_res)
    s = substitution_matrix(ref_res, query_res)

    H, T = score_matrix(s+d)
    alignment = traceback(H, T)
    return np.max(H), alignment

def protalign(ref, query, min_score=0):
    """
    Assign protein residues positioned similarly in sequence and 3D space
    the same chain and residue number.
    """
    mapping = {}
    for query_chain in query.chain:
        best_score, best_alignment, best_chain = min_score, None, None
        for ref_chain in ref.chain:
            score, alignment = align_chains(ref_chain, query_chain)

            if score > best_score:
                best_score = score
                best_alignment = alignment
                best_chain = ref_chain
        
        if best_alignment is not None:
            ref_res = get_residues(ref_chain, centroid=False)
            query_res = get_residues(query_chain, centroid=False)
            for i, j in zip(*best_alignment):
                if i is not None and j is not None:
                    mapping[query_res[j-1]] = ref_res[i-1]
    return mapping

def get_mapping(ref, query):
    prot_mapping = protalign(ref, query)
    het_mapping = hetalign(ref, query)
    return prot_mapping
