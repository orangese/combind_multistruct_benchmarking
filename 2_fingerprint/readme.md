Goal: generate a fingerprint for every pose.

## 1) Fingerprint the crystallized pose

`./crystal_main.py <receptor> <output_dir>`

`output_dir` will show up in `/scratch/PI/rondror/docking_data/<receptor>/`.

For example:

`./crystal_main.py TRMD crystal_ifps`

## 2) Fingerprint the docked poses

`./glides_main.py <receptor> <output_dir>`

