
echo Ligand with fused rings
$SCHRODINGER/run ~/combind/ifp/fp.py -poses 1 -mode pv -input_file inputs/2HVC_lig-to-1T5Z_pv.maegz -output_file temp
grep '^5' temp > outputs/2HVC_lig-to-1T5Z.fp
rm temp


echo Ligand with fused ring and 2 disjoint rings
$SCHRODINGER/run ~/combind/ifp/fp.py -poses 1 -mode pv -input_file inputs/4AQC_lig-to-2B7A_pv.maegz -output_file temp
grep '^5' temp > outputs/4AQC_lig-to-2B7A.fp
rm temp
