
echo case 1
$SCHRODINGER/run ~/combind/mcss/mcss.py INIT CHEMBL3092760_lig CHEMBL3310115_lig inputs/CHEMBL3092760_lig.mae inputs/CHEMBL3310115_lig.mae outputs/CHEMBL3092760_lig-CHEMBL3310115_lig.csv ~/combind/mcss/custom_types/mcss14.typ
diff outputs/CHEMBL3092760_lig-CHEMBL3310115_lig.csv expected/CHEMBL3092760_lig-CHEMBL3310115_lig.csv

echo case 2
$SCHRODINGER/run ~/combind/mcss/mcss.py INIT 3AUQ_lig CHEMBL3264162_lig inputs/3AUQ_lig.mae inputs/CHEMBL3264162_lig.mae outputs/3AUQ_lig-CHEMBL3264162_lig.csv ~/combind/mcss/custom_types/mcss14.typ
diff outputs/3AUQ_lig-CHEMBL3264162_lig.csv expected/3AUQ_lig-CHEMBL3264162_lig.csv

echo case 3
$SCHRODINGER/run ~/combind/mcss/mcss.py INIT CHEMBL210836_lig CHEMBL424880_lig inputs/CHEMBL210836_lig.mae inputs/CHEMBL424880_lig.mae outputs/CHEMBL210836_lig-CHEMBL424880_lig.csv ~/combind/mcss/custom_types/mcss14.typ
diff outputs/CHEMBL210836_lig-CHEMBL424880_lig.csv expected/CHEMBL210836_lig-CHEMBL424880_lig.csv
