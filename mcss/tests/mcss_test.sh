
echo MCSS INIT with atomic weight specified
echo TODO
echo MCSS RMSD with atomic weight specified
echo TODO

echo No matches
$SCHRODINGER/run ~/combind/2_ifp/mcss_main.py RMSD CHEMBL3092760_lig CHEMBL3310115_lig inputs/CHEMBL3092760_lig-to-4IB4_pv.maegz inputs/CHEMBL3310115_lig-to-4IB4_pv.maegz outputs/CHEMBL3092760_lig-CHEMBL3310115_lig-4IB4-glide12.csv 100  "CHEMBL3092760_lig,CHEMBL3310115_lig,19,21,19,18,ccCcc(cccc)ccccC(=O)[N+]=C(N)N,ccCcc(cccc)ccccC(=O)[N+]=C(N)N"

$SCHRODINGER/run ~/combind/2_ifp/mcss_main.py RMSD CHEMBL3092760_lig CHEMBL3310115_lig inputs/CHEMBL3092760_lig-to-4IB4_pv.maegz inputs/CHEMBL3310115_lig-to-4IB4_pv.maegz outputs/CHEMBL3092760_lig-CHEMBL3310115_lig-4IB4-glide12.csv 100  "CHEMBL3092760_lig,CHEMBL3310115_lig,19,21,19,18,ccCcc(cccc)ccccC(=O)[N]=C(N)N,ccCcc(cccc)ccccC(=O)[N+]=C(N)N"

echo MCSS RMSD with out of order smiles expression
$SCHRODINGER/run ~/combind/2_ifp/mcss_main.py RMSD 3AUQ_lig CHEMBL3264162_lig inputs/3AUQ_lig-to-1DB1_pv.maegz inputs/CHEMBL3264162_lig-to-1DB1_pv.maegz outputs/3AUQ_lig-CHEMBL3264162_lig-1DB1-glide12.csv 10  "3AUQ_lig,CHEMBL3264162_lig,30,33,28,30,CCCCCC(C)C1CCC(C12C)C(=CCC2)C#CC3=CC(O)C(C)C(C3)O,CCCCCC(C)C1CCC(C12C)C(CCC2)=CC=C3CC(O)C(=C)C(O)C3"
diff expected/3AUQ_lig-CHEMBL3264162_lig-1DB1-glide12.csv outputs/3AUQ_lig-CHEMBL3264162_lig-1DB1-glide12.csv

echo 'MCSS RMSD does not write if no matches (expect no such file message)'
$SCHRODINGER/run ~/combind/2_ifp/mcss_main.py RMSD CHEMBL210836_lig CHEMBL424880_lig inputs/CHEMBL210836_lig-to-1BZC_pv.maegz inputs/CHEMBL424880_lig-to-1BZC_pv.maegz outputs/CHEMBL210836_lig-CHEMBL424880_lig-1BZC-glide12.csv 10  "CHEMBL210836_lig,CHEMBL424880_lig,35,36,35,39,c1ccccc1S(=O)(=O)NC(c(n2)nc(c23)cccc3)Cc4ccc(cc4)C(S5(=O)=O)CC(=N5)[O-],c1ccccc1S(=O)(=O)NC(c(n2)nc(c23)cccc3)Cc4ccc(cc4)C(S5(=O)=O)CC(=N5)[O-]"
ls outputs/outputs/CHEMBL210836_lig-CHEMBL424880_lig-1BZC-glide12.csv
