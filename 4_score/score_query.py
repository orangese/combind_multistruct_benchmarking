import os
import sys
import numpy

script_path, q, s = sys.argv
code_path = '/'.join(script_path.split('/')[:-2])

for i in ['1_dock','2_fp','3_analyze','4_score']:
    sys.path.append(code_path+'/'+i)

from sc_containers import ScoreContainer

full_output_dir = os.getcwd().split('/')
data_dir = '/'.join(full_output_dir[:7])
p = full_output_dir[7]
output_dir = '/'.join(full_output_dir[8:])

sc = ScoreContainer(data_dir, output_dir, p, s)
sc.compute_results(q)


