import matplotlib.pyplot as plt
import matplotlib.cm as cm
from sklearn.cluster import KMeans
import numpy as np

class CrystalCluster:
    def __init__(self, crystals, n, i_indices):
        self.crystals = crystals
        self.i_indices = i_indices

        alpha_sort = sorted(crystals.keys())
        self.kmeans = KMeans(n_clusters=n, random_state=0).fit(self.get_fp_matrix(alpha_sort)[0])
        
        self.sorted_ligs = sorted(alpha_sort, key=lambda x:self.kmeans.labels_[alpha_sort.index(x)])
        
        self.lig_clusters = { 
            i : [l for l in alpha_sort if self.kmeans.labels_[alpha_sort.index(l)] == i] for i in range(n) }
        
        self.fp_mat, self.all_i = self.get_fp_matrix(self.sorted_ligs)

    def show_clusters_together(self, num_i=15):
        self.interaction_heatmap(self.fp_mat[:,:num_i], self.sorted_ligs, self.all_i[:num_i], self.kmeans.labels_)
        
    def show_cluster_centers(self, num_i=15):
        self.interaction_heatmap(self.kmeans.cluster_centers_[:,:num_i], 
                                 [i for i in range(max(self.kmeans.labels_) + 1)], self.all_i[:num_i])
        
    def show_clusters_individually(self, num_i=15):
        for i, c_ligs in self.lig_clusters.items():
            fp_mat, all_i = self.get_fp_matrix(c_ligs)
            self.interaction_heatmap(fp_mat[:,:num_i], c_ligs, all_i[:num_i])
        
    def get_fp_vectors(self, ligs):
        all_interactions = {}
        fp_map = {}
        for s in ligs:
            if s not in fp_map:
                fp_map[s] = {}
            for r in self.crystals[s]:
                for i in self.i_indices:
                    if self.crystals[s][r][i] != 0:
                        fp_map[s][(r,i)] = self.crystals[s][r][i]
                        all_interactions[(r,i)] = all_interactions.get( (r,i), 0) + fp_map[s][(r,i)]

        all_interactions = sorted(all_interactions.keys(), key=lambda x:-all_interactions[x])
        return all_interactions, {s: [fp_map[s].get(x, 0) for x in all_interactions] for s in ligs}
        
    def get_fp_matrix(self, sorted_ligs):
        all_i, fp_vectors = self.get_fp_vectors(sorted_ligs)
        fp_matrix = np.zeros(( len(sorted_ligs), len(fp_vectors[sorted_ligs[0]]) ))
        for i, l in enumerate(sorted_ligs):
            fp_matrix[i,:] = fp_vectors[l][:]
        return fp_matrix, all_i
        
    def interaction_heatmap(self, A, structs, res, cluster_labels=None):
        fig, ax = plt.subplots()
        
        sq_size = 0.3
        fig.set_size_inches(sq_size*len(res), sq_size*len(structs), forward=True)
        
        def i_matrix(A, res, i):
            aa = np.zeros(A.shape)
            for j, r in enumerate(res):
                if r[1] == i: aa[:,j] = A[:,j]
                else: aa[:,j] = np.nan*A[:,j]
            return aa
        
        colors = { 0:cm.Purples, 1:cm.Blues, 2:cm.Greens, 3:cm.Reds, 4:cm.Greys }
        for j,i in enumerate(self.i_indices):        
            ax.matshow(i_matrix(A, res, i), cmap=colors[j%len(colors.keys())], vmin=0, vmax=np.max(A))

        # put the major ticks at the middle of each cell
        ax.set_xticks(np.arange(A.shape[1]), minor=False)
        ax.set_yticks(np.arange(A.shape[0]), minor=False)
        ax.xaxis.tick_top()

        if cluster_labels is not None:
            for c in range(max(cluster_labels)):
                num_in_c = sum([1 for i in cluster_labels if i <= c]) - 0.5
                ax.plot(ax.get_xlim(), [num_in_c, num_in_c], linewidth = 4, c="k")

        ax.set_xticklabels(res, minor=False, rotation = 'vertical')
        ax.set_yticklabels(structs, minor=False)
        plt.xlabel('interactions')
        plt.ylabel('crystal structures')

        plt.show()