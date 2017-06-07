############# Glide Cumulative Based prediction

#For a given ligand -> grid, take the potential poses and rerank them by:
#For every other ligand2 -> grid, find the best overlap and add it to a cluster
#Reweight the potential poses by the sum of each cluster

# auxiliary funciton to make it work
def glide_predict_helper(args):
    return glide_predict(*args)

def glide_predict(i, j, ligand, grid):
    clusters = []
    #Make a cluster of best TC FP overlaped poses for every pose generated from docking the given ligand -> grid
    for cluster in map(Cluster, glides[ligand][grid].poses.values()):
        for ligand2 in (all_gridstructs + ligstructs):
            if ligand2 == ligand: continue
            best = (0, None)
            for pose in glides[ligand2][grid].poses.values(): # compare poses of ligand2 to the given pose of ligand and find the best-overlapping one
                score = pose.fp.jaccard_index(cluster.lead.fp)
                if score > best[0] or best[1] is None:
                    best = (score, pose)
            cluster.add_pose(best[1])
        clusters += [cluster]
    return i, j, sorted(clusters, key = lambda x: x.jaccard_index()) #Sort our doced ligand -> grid poses by total cumulative overlap

def run_glides(num_poses):
    pool = Pool(ncpus)

    A = np.zeros((len(gridstructs), len(all_gridstructs)))

    threadingArguments = [(tgrid[0], tall[0], tgrid[1], tall[1]) for tgrid, tall in itertools.product(enumerate(gridstructs), enumerate(all_gridstructs))]

    for result in pool.imap_unordered(glide_predict_helper, threadingArguments):
        i, j, data = result
        A[i,j] = data[-1].lead.rmsd

    print 'Glide comparison: max of A is: ' + str(np.amax(A))
    heatmap(A, 'jaccard_glide_'+str(num_poses))
    fnm = make_filename('jaccard_glide_'+str(num_poses)+'.npy')
    np.save(fnm, A)
    return A


# In[ ]:

#First, run glide with all 150 poses
print('Running 150 poses')
                                  
# In[ ]:

#Cull the number of poses that we use to predict
for num_pose in [100, 50, 25]:
    for tKey in glides.keys():
        for tKey2 in glides[tKey].keys():
            numKeys = max(glides[tKey][tKey2].poses.keys())
            maxKey = min(num_pose, numKeys)
            glides[tKey][tKey2].poses = {your_key: glides[tKey][tKey2].poses[your_key] for your_key in range(maxKey)}
    print('Running ' + str(num_pose) + ' poses')
    glide = run_glides(num_pose)


    def __str__(self):
        return str(self.overlap())

