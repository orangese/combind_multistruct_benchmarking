Goal: generate a fingerprint for every pose.

## 1) Fingerprint the crystallized pose

`./crystal_main.py <output_dir> <receptor 1> <receptor 2> ...`

`output_dir` will show up in `/scratch/PI/rondror/docking_data/<receptor>/ifp`.

For example:

`./crystal_main.py crystal_ifps B1AR B2AR TRMD`

## 2) Fingerprint the docked poses

`./glide_main.py <output_dir> <receptor_1> <receptor_2> ...`

