
# coding: utf-8

# In[4]:

from multiprocessing import Pool
import os
import sys
import tqdm
import itertools
from math import exp
from tqdm import tnrange, tqdm_notebook, tqdm


# In[5]:

ncpus = 4 #DEFINE -c CPUS_PER_NODE


# In[6]:

sys.path.append('/scratch/PI/rondror/docking/workflow_thomas/fingerprint')


# In[7]:

os.chdir("/scratch/PI/rondror/docking/data/TGFB1")


# In[8]:

from fingerprint import FuzzyFingerPrint
from rmsd import RMSD
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import cm as CM
import numpy as np
import sys
import os.path


# In[9]:

#Container for information about a Pose
class Pose:
    def __init__(self, rmsd, fp, num, gscore):
        self.rmsd = rmsd
        self.fp = fp
        self.num = num
        self.gscore = gscore

    def __str__(self):
        return ';'.join(map(str, [self.num, self.rmsd, self.fp, self.gscore]))

#Container for information about a ligand
class Ligand:
    def __init__(self, crystal):
        self.crystal = crystal
        self.poses = {}

    def add_pose(self, pose, pose_num):
        self.poses[pose_num] = pose

    def __str__(self):
        return str(self.crystal)

#Returns the location of where we should output our NumPy arrays to one we have performed the analysis
def make_filename(fnm):
  if ONEGRID:
    fnm = gridstructs[0] + '_' + fnm
  if output_dir != "":
    fnm = output_dir + '/' + fnm
  return fnm


# In[10]:

sys.argv = 'file.py TGFB1 all ligands grids glide rmsd.csv crystalfifps/fuzzy_ifp.fp newfifp sar_analysis_thomas'.split(" ")


# In[11]:

#Read in command line inputs
dataDirectory = "/scratch/PI/rondror/docking/data/"

data_set_dir = sys.argv[1]
sar_structures = sys.argv[2]
ligands_dir = sys.argv[3]
grids_dir = sys.argv[4]
glide_dir = sys.argv[5]
rmsd_file = sys.argv[6]
crystal_fp_file = sys.argv[7]
docking_fp_dir = sys.argv[8]
output_dir = sys.argv[9]

#Make sure things were inputted correctly
if data_set_dir[-1] != '/':
  data_set_dir += '/'

if ligands_dir[-1] != '/':
  ligands_dir += '/'

if grids_dir[-1] != '/':
  grids_dir += '/'

if glide_dir[-1] != '/':
  glide_dir += '/'

if '.csv' not in rmsd_file:
  rmsd_file += '.csv'

if docking_fp_dir[-1] != '/':
  docking_fp_dir += '/'

if output_dir[-1] != '/':
  output_dir += '/'

#All directories should be full paths
data_set_dir = dataDirectory + data_set_dir
ligands_dir = data_set_dir + ligands_dir
grids_dir = data_set_dir + grids_dir
glide_dir = data_set_dir + glide_dir
docking_fp_dir = data_set_dir + docking_fp_dir
output_dir = data_set_dir + output_dir


# In[12]:

#Get a list of all grid structures that we have
all_gridstructs = [d.lower() for d in os.listdir(grids_dir) if os.path.isdir(os.path.join(grids_dir, d))]

if '5e8w' in all_gridstructs:
    all_gridstructs.remove('5e8w')


# In[13]:

#Get a list of all ligands that we have
ligstructs = [d for d in os.listdir(ligands_dir) if os.path.isfile(os.path.join(ligands_dir, d))]
ligstructs = map(lambda x: x.split("_")[0], ligstructs)

if '5e8w' in ligstructs:
    ligstructs.remove('5e8w')


# In[14]:

#gridstructs contains a list of structure that we want to perform SAR analysis on
ONEGRID = False
if sar_structures.strip().lower() == 'all':
    gridstructs = all_gridstructs
else:
    gridstructs = [sys.argv[1]]
    ONEGRID = True

if '3qak' in all_gridstructs:
    all_gridstructs.remove('3qak')
    
#Sort our ligands by the order in which we perform SAR analysis
all_gridstructs_upper = map(lambda x : x.upper(), all_gridstructs)
ligstructs.sort(key=lambda lig: all_gridstructs_upper.index(lig.upper()) if lig.upper() in all_gridstructs_upper else len(all_gridstructs_upper) + 1)


# In[15]:

# Ligand -> Grid Structure -> List of Glide Scores sorted by pose number (Lowest Scores to Highest Scores)
gscores = {}

for ligand in all_gridstructs:
    gscores[ligand] = {}
    for grid in all_gridstructs:
        gscores[ligand][grid] = []
        
        #Go through all possible permutations for file capitalization
        #NOTE: Files still have to be in the {}_ligand-to-{} format!
        cap_permutations = [(ligand, grid), (ligand.upper(), grid.upper()), (ligand.upper(), grid.lower()), (ligand.lower(), grid.upper()), (ligand.lower(), grid.lower())]
        for tLigand, tGrid in cap_permutations:
            fnm = glide_dir + "{}_ligand-to-{}/{}_ligand-to-{}.rept".format(tLigand, tGrid, tLigand, tGrid)
            try:
                s_file = open(fnm)
                break
            except:
                pass
        
        line = s_file.readline().strip().split()
        while not line or line[0] != '1' or (len(line) != 19 and (len(line) > 1 and line[1] != "1" and len(line) != 18)): # hack - why is title blank sometimes?
            line = s_file.readline().strip().split()
        while line:
            # Rank', 'Title', 'Lig#', 'Score', 'GScore'
            rank, title, lig, score, gscore = line[:5]
            gscores[ligand][grid] += [float(gscore)]
            line = s_file.readline().strip().split()


# In[16]:

# Ligand -> Grid Structure -> List of RMSDs ordered by GLIDE Pose Number
rmsds = {}
#for line in open('rmsd_table.csv'):
for line in open('rmsd.csv'):
    n, data = line.strip().split(':')
    ligand, grid = n.split('-to-')
    ligand = ligand.split('_')[0]
    grid = grid.split('_')[0]
        
    if grid.lower() not in rmsds: rmsds[grid.lower()] = {}
    rmsds[grid.lower()][ligand.lower()] = RMSD.read(data)


# In[17]:

glides = {}
for ligand in all_gridstructs:
    glides[ligand] = {}
    for grid in all_gridstructs:
        glides[ligand][grid] = Ligand(None)
        
        #Go through all possible permutations for file capitalization
        #NOTE: Files still have to be in the {}_ligand-to-{} format!
        fnm = docking_fp_dir + "{}_ligand-to-{}.fp".format(ligand, grid)
        try:
            s_file = open(fnm)
        except:
            pass

        for pose_num, line in enumerate(s_file):
            try:
                ifp = line.strip()
            except Exception as e:
                print("Error: No fingerprint generated for pose " + line.strip() + " in " + fnm)
                break
            
            if pose_num < 50: #ONLY IMPORT 50 POSES, REMOVE LATER
                glides[ligand][grid].add_pose(Pose(rmsds[grid][ligand].get_rmsd(pose_num),
                                                   FuzzyFingerPrint.compact_parser(ifp, ligand), pose_num,
                                                   gscores[ligand][grid][pose_num]), pose_num)


# In[15]:

############### Shared Classes

### Compare RMSDs to Fingerprint scores

def rmsd_to_fp(ex):
    X, Y = [], []
    for pose in glides[ex][ex].poses:
        pose = glides[ex][ex].poses[pose]
        print pose.num, pose.rmsd, pose.fp.overlap(crystals[ex].fp)
        X += [pose.rmsd]
        Y += [pose.fp.overlap(crystals[ex].fp)]
    plt.scatter(X, Y)
    plt.savefig('test.png')
    
class Cluster:
    def __init__(self, poses): #poses is a list of Pose objects
        self.floatsPerResidue = 6 #Could change if we change the fingerprint
        self.poses = poses
        self.residues = list(set().union(*[pose.fp.feats.keys() for pose in self.poses])) #Cache all of the residues
        self.numpyFPs = []
        self.score = None
        
        #Setup the numpy vectorized FPs
        for pose in self.poses:
            lengthFPVector = self.floatsPerResidue * len(self.residues)
            poseFP = np.zeros(lengthFPVector)
            rawPoseFP = pose.fp.feats
            
            for residueIndex, residue in enumerate(self.residues):
                if residue in rawPoseFP:
                    poseFP[residueIndex * self.floatsPerResidue:(residueIndex+1) * self.floatsPerResidue] = rawPoseFP[residue]
                    
            self.numpyFPs.append(poseFP)

    def get_individual_overlap(self, a):
        if np.sign(a[0]) == np.sign(a[1]): #If they are the same sign
            return np.minimum(np.abs(a[0]), np.abs(a[1]))
        return 0
    
    def overlap(self):
        totalOverlap = 0
        
        for indOne, fpOne in enumerate(self.numpyFPs):
            for indTwo, fpTwo in enumerate(self.numpyFPs):
                if indOne != indTwo:
                    newFP = np.vstack((fpOne, fpTwo))
                    totalOverlap += np.sum(np.apply_along_axis(self.get_individual_overlap, 0, newFP))
                
        return totalOverlap
    
    def overlap_per_pair(self):
        totalOverlap_ind = []
        ind_track = []
        
        for indOne, fpOne in enumerate(self.numpyFPs):
            for indTwo, fpTwo in enumerate(self.numpyFPs):
                if indOne != indTwo:
                    newFP = np.vstack((fpOne, fpTwo))
                    totalOverlap_ind.append(np.sum(np.apply_along_axis(self.get_individual_overlap, 0, newFP)))
                    ind_track.append((indOne, indTwo))
                
        return totalOverlap_ind, ind_track

    def __str__(self):
        return str(self.overlap())

def heatmap(A, plotTitle):
    fig, ax = plt.subplots()
    heatmap = ax.pcolor(A, cmap=CM.jet, vmin=0, vmax=10.25)
    plt.colorbar(heatmap)

    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(A.shape[1]) + 0.5, minor=False)
    ax.set_yticks(np.arange(A.shape[0]) + 0.5, minor=False)
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    #labels
    column_labels = all_gridstructs
    row_labels = all_gridstructs
    ax.set_xticklabels(column_labels, minor=False, rotation = 'vertical')
    ax.set_yticklabels(row_labels, minor=False)
    ax.plot(ax.get_xlim(), ax.get_ylim()[::-1], linewidth = 4, c="m")
    plt.title(plotTitle)
    plt.savefig(plotTitle + '.png')
    plt.show()


# In[16]:

################ Top pose

def crystal_energy(exclude = None):
    out = {}
    for crystal in crystals:
        if crystal == exclude: continue
        fp = crystals[crystal].fp
        for resi in fp.residues():
            common = fp.to_common(resi)
            if not common in out: out[common] = [0] * 7 # could use defaultdict instead? This vector represents the interactions that a particular residue can make: T-stack, Pi=pi, cat-pi, salt, hydro, donor, accept (see PPT slide 22)
            out[common] = [i+j for i, j in zip(out[common], fp.entry(resi))]

    return out

# Top scoring GLIDE pose
def top_pose():
    A = np.zeros((len(gridstructs), len(all_gridstructs)))
    for i, grid in enumerate(gridstructs):
        for j, ligand in enumerate(all_gridstructs):
            A[i, j] = rmsds[grid][ligand].get_rmsd(0)
    #heatmap(A, 'top_rmsd.png')
    return A

# Heatmap of the tanimoto coefficients, instead of RMSDs
# why does it look weird?
def top_pose_tc():
    A = np.zeros((len(gridstructs), len(all_gridstructs)))
    for i, grid in enumerate(gridstructs):
        for j, ligand in enumerate(all_gridstructs):
            A[i, j] = crystals[ligand].fp.tanimoto_coef(glides[ligand][grid].poses[0].fp)
    #heatmap(A, 'top_tc.png')
    return A

# Best out of all generated poses
def best_pose():
    A = np.zeros((len(gridstructs), len(all_gridstructs)))
    for i, grid in enumerate(gridstructs):
        for j, ligand in enumerate(all_gridstructs):
            A[i, j] = rmsds[grid][ligand].best()
    print 'Max of best: ' + str(np.amax(A))
    return A

# Best out of all generated poses
def best_pose_index():
    A = np.zeros((len(gridstructs), len(all_gridstructs)))
    for i, grid in enumerate(gridstructs):
        for j, ligand in enumerate(all_gridstructs):
            A[i, j] = rmsds[grid][ligand].getMinPose()
    print 'Max of best: ' + str(np.amax(A))
    return A

# again, tanimoto coefficients...
def best_tc():
    A = np.zeros((len(gridstructs), len(all_gridstructs)))
    for i, grid in enumerate(gridstructs):
        for j, ligand in enumerate(all_gridstructs):
            A[i, j] = max(pose.fp.tanimoto_coef(crystals[ligand].fp) for pose in glides[ligand][grid].poses.values())
    #heatmap(A, 'best_tc.png')
    fnm = make_filename('best_tc.npy')
    np.save(fnm, A)
    return A


# In[17]:

from tqdm import tqdm
from random import random
from copy import copy


# In[ ]:

from random import randint

#Initialize Variables
old_rmsds = []
new_rmsds = []
labels = []

for grid in tqdm(all_gridstructs):
    if grid == '2vt4':
        continue
    
    best_poses = []
    best_poses_index = [0] * len(all_gridstructs) #Start off with [0,...,0] as our intial set
    best_overlap = None
    overlap_history = []
    pose_history = []

    #Initialize the starting poses to just be the poses at index 0
    for index, ligand in enumerate(all_gridstructs):
        best_poses.append(glides[ligand][grid].poses[best_poses_index[index]])

    cluster = Cluster(best_poses)
    best_overlap = cluster.overlap()

    #Start MCMC Procedure
    num_trial = 0
    for num_trial in tqdm(range(100)):
        #Find the maximum overlap for randomIndex's ligand
        randomIndex = randint(0,len(best_poses)-1) 
        temp_poses = [glides[ligand][grid].poses[best_poses_index[index]]
                      for index, ligand in enumerate(all_gridstructs)]
        best_temp_poses = None
        best_temp_overlap = -1
        best_random_index = -1

        randomIndexPoses = glides[all_gridstructs[randomIndex]][grid].poses
        for index in range(len(randomIndexPoses)):
            temp_poses[randomIndex] = randomIndexPoses[index]
            cluster = Cluster(temp_poses)
            cluster_overlap = cluster.overlap()

            #Find the best overlap that is not the original overlap that we started with
            if cluster_overlap > best_temp_overlap and cluster_overlap != best_overlap:
                best_temp_poses = temp_poses
                best_temp_overlap = cluster_overlap
                best_random_index = index

        #If the best overlap is lower, change it - if not, still change it with random probability
        print("Best Overlap: " + str(best_overlap))
        print("Temp Overlap: " + str(best_temp_overlap))
        if best_temp_overlap > best_overlap or (np.exp(float(best_temp_overlap - best_overlap)/200) >= random()):
            best_poses = best_temp_poses
            best_poses_index[randomIndex] = best_random_index
            best_overlap = best_temp_overlap

        pose_history.append((copy(best_poses_index), best_overlap))
        overlap_history.append(best_overlap)
        
    #After MCMC is finished, append the best index (scored by overlap) so that we can plot later
    best_indexes = sorted(pose_history, key=lambda x:x[1], reverse=True)[0][0]
        
    for index, temp_index in enumerate(best_indexes):
        new_rmsds.append(rmsds[grid][all_gridstructs[index]].data[temp_index])

    for index, temp_index in enumerate([0] * len(all_gridstructs)):
        old_rmsds.append(rmsds[grid][all_gridstructs[index]].data[temp_index])
    
    for index, temp_index in enumerate([0] * len(all_gridstructs)):
        labels.append(str(all_gridstructs[index]) + "_ligand-to-" + grid +"_grid")


# In[ ]:

import pickle
 
with open('new_rmsds.pkl', 'wb') as f:
    pickle.dump(new_rmsds, f)
    
with open('old_rmsds.pkl', 'wb') as f:
    pickle.dump(old_rmsds, f)
    
with open('labels.pkl', 'wb') as f:
    pickle.dump(labels, f)


# In[ ]:

import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')

plt.plot([0, 12], [0,12])
plt.plot(new_rmsds, old_rmsds, ".")
plt.xlabel("MCMC RMSD")
plt.ylabel("Glide RMSD")
plt.title("CDK2 Glide vs. MCMC RMSD")
plt.savefig('glide_vs_mcmc.png')
plt.show()


# In[ ]:

import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')

texts = []
plt.figure(figsize=(20,15))
plt.plot([0, 12], [0,12])
plt.plot(new_rmsds, old_rmsds, ".")
for num, new_old in enumerate(zip(new_rmsds, old_rmsds)):
    new = new_old[0]
    old = new_old[1]
    texts.append(plt.text(new, old, labels[num]))
    
plt.xlabel("MCMC RMSD")
plt.ylabel("Glide RMSD")
plt.title("CDK2 Glide vs. MCMC RMSD")
plt.savefig('glide_vs_mcmc_labeled.png')
plt.show()


# In[ ]:

bestRMSD = []
topRMSD = []

for grid in all_gridstructs:
    for ligand in all_gridstructs:
        rmsdList = rmsds[grid][ligand].data
        bestRMSD.append(min(rmsdList[0:50]))
        topRMSD.append(rmsdList[0])


# In[ ]:

import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')

plt.plot([0, 12], [0,12])
plt.plot(bestRMSD, topRMSD, ".")
plt.xlabel("Best RMSD")
plt.ylabel("Top RMSD")
plt.title("CDK2 Top vs. Best RMSD")
plt.savefig('top_vs_best.png')
plt.show()

