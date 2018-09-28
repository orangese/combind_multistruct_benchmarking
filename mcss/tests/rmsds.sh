
echo MCSS RMSD with atomic weight specified
echo TODO

echo MCSS size decreases mcss14 fails
$SCHRODINGER/run ~/combind/mcss/mcss.py RMSD CHEMBL186288_lig CHEMBL187750_lig inputs/CHEMBL186288_lig-to-1AQ1_pv.maegz inputs/CHEMBL187750_lig-to-1AQ1_pv.maegz outputs/CHEMBL186288_lig-CHEMBL187750_lig.init.csv ~/combind/mcss/custom_types/mcss14.typ outputs/CHEMBL186288_lig-CHEMBL187750_lig-1AQ1-glide12.csv 10  "CHEMBL186288_lig,CHEMBL187750_lig,27,27,23,26,c1ccccc1Nc2nccc(n2)-c(cn3)c(n34)ccc(O)n4,c1ccccc1Nc2nccc(n2)-c(cn3)c(n34)ccc([O-])n4"

echo MCSS size decreases mcss15
$SCHRODINGER/run ~/combind/mcss/mcss.py RMSD CHEMBL186288_lig CHEMBL187750_lig inputs/CHEMBL186288_lig-to-1AQ1_pv.maegz inputs/CHEMBL187750_lig-to-1AQ1_pv.maegz outputs/CHEMBL186288_lig-CHEMBL187750_lig.init.csv ~/combind/mcss/custom_types/mcss15.typ outputs/CHEMBL186288_lig-CHEMBL187750_lig-1AQ1-glide12.csv 10  "CHEMBL186288_lig,CHEMBL187750_lig,27,27,23,26,c1ccccc1Nc2nccc(n2)-c(cn3)c(n34)ccc(O)n4,c1ccccc1Nc2nccc(n2)-c(cn3)c(n34)ccc([O-])n4"

echo Mysteriously does not work with old schrodinger.structutils.analyze.evaluate_smarts
$SCHRODINGER/run ~/combind/mcss/mcss.py RMSD 3PSD_lig CHEMBL1079785_lig inputs/3PSD_lig-to-1UWH_pv.maegz inputs/CHEMBL1079785_lig-to-1UWH_pv.maegz outputs/3PSD_lig-CHEMBL1079785_lig.init.csv ~/combind/mcss/custom_types/mcss14.typ outputs/3PSD_lig-CHEMBL1079785_lig-1UWH-glide12.csv 10  "3PSD_lig,CHEMBL1079785_lig,29,37,21,23,c1cnccc1-c2nn(CCC[N+])cc2-c3ccccc3,c1cnccc1-c(n2)c(cn2cccn)-c3ccccc3"

echo No matches
$SCHRODINGER/run ~/combind/mcss/mcss.py RMSD CHEMBL3092760_lig CHEMBL3310115_lig inputs/CHEMBL3092760_lig-to-4IB4_pv.maegz inputs/CHEMBL3310115_lig-to-4IB4_pv.maegz outputs/CHEMBL3092760_lig-CHEMBL3310115_lig.csv ~/combind/mcss/custom_types/mcss14.typ outputs/CHEMBL3092760_lig-CHEMBL3310115_lig-4IB4-glide12.csv 10  "CHEMBL3092760_lig,CHEMBL3310115_lig,19,21,19,18,ccCcc(cccc)ccccC(=O)[N+]=C(N)N,ccCcc(cccc)ccccC(=O)[N+]=C(N)N"
diff expected/CHEMBL3092760_lig-CHEMBL3310115_lig-4IB4-glide12.csv outputs/CHEMBL3092760_lig-CHEMBL3310115_lig-4IB4-glide12.csv

echo MCSS RMSD with out of order smiles expression
$SCHRODINGER/run ~/combind/mcss/mcss.py RMSD 3AUQ_lig CHEMBL3264162_lig inputs/3AUQ_lig-to-1DB1_pv.maegz inputs/CHEMBL3264162_lig-to-1DB1_pv.maegz outputs/3AUQ_lig-CHEMBL3264162_lig.csv ~/combind/mcss/custom_types/mcss14.typ outputs/3AUQ_lig-CHEMBL3264162_lig-1DB1-glide12.csv 10 "3AUQ_lig,CHEMBL3264162_lig,30,33,28,30,CCCCCC(C)C1CCC(C12C)C(=CCC2)C#CC3=CC(O)C(C)C(C3)O,CCCCCC(C)C1CCC(C12C)C(CCC2)=CC=C3CC(O)C(=C)C(O)C3"
diff expected/3AUQ_lig-CHEMBL3264162_lig-1DB1-glide12.csv outputs/3AUQ_lig-CHEMBL3264162_lig-1DB1-glide12.csv

echo Tautomers differ between prepared ligands and glide poses
$SCHRODINGER/run ~/combind/mcss/mcss.py RMSD CHEMBL210836_lig CHEMBL424880_lig inputs/CHEMBL210836_lig-to-1BZC_pv.maegz inputs/CHEMBL424880_lig-to-1BZC_pv.maegz outputs/CHEMBL210836_lig-CHEMBL424880_lig.csv ~/combind/mcss/custom_types/mcss14.typ outputs/CHEMBL210836_lig-CHEMBL424880_lig-1BZC-glide12.csv 10  "CHEMBL210836_lig,CHEMBL424880_lig,35,36,35,39,c1ccccc1S(=O)(=O)NC(c(n2)nc(c23)cccc3)Cc4ccc(cc4)C(S5(=O)=O)CC(=N5)[O-],c1ccccc1S(=O)(=O)NC(c(n2)nc(c23)cccc3)Cc4ccc(cc4)C(S5(=O)=O)CC(=N5)[O-]"
diff expected/CHEMBL210836_lig-CHEMBL424880_lig-1BZC-glide12.csv outputs/CHEMBL210836_lig-CHEMBL424880_lig-1BZC-glide12.csv


echo Ketone Oxygen to oxyanion conversion mcss14
$SCHRODINGER/run ~/combind/mcss/mcss.py RMSD CHEMBL2385551_lig CHEMBL2385552_lig inputs/CHEMBL2385551_lig-to-1KV1_pv.maegz inputs/CHEMBL2385552_lig-to-1KV1_pv.maegz outputs/CHEMBL2385551_lig-CHEMBL2385552_lig.init.csv ~/combind/mcss/custom_types/mcss14.typ outputs/CHEMBL2385551_lig-CHEMBL2385552_lig-1KV1-glide12.csv 10  "CHEMBL2385551_lig,CHEMBL2385552_lig,27,27,24,26,CSc1ccc(cc1)-c2nc(-c(cc)cc)c(n2)-c3ccccc3,C[S+]c1ccc(cc1)-c2nc(C(=CC)CC)c(n2)-c3ccccc3"

echo Ketone Oxygen to oxyanion conversion mcss15
$SCHRODINGER/run ~/combind/mcss/mcss.py RMSD CHEMBL2385551_lig CHEMBL2385552_lig inputs/CHEMBL2385551_lig-to-1KV1_pv.maegz inputs/CHEMBL2385552_lig-to-1KV1_pv.maegz outputs/CHEMBL2385551_lig-CHEMBL2385552_lig.init.csv ~/combind/mcss/custom_types/mcss15.typ outputs/CHEMBL2385551_lig-CHEMBL2385552_lig-1KV1-glide12.csv 10  "CHEMBL2385551_lig,CHEMBL2385552_lig,27,27,25,0,CS(=O)c1ccc(cc1)-c2nc(-c(cc)cc)c(n2)-c3ccccc3,C[S+]([O-])c1ccc(cc1)-c2nc(C(CC)=CC)c(n2)-c3ccccc3"
