
echo 'Ligand with imine group'
$SCHRODINGER/run ~/combind/ifp/fp.py -poses 1 -mode pv -input_file inputs/1F5L_lig-to-1C5X_pv.maegz -output_file inputs/1F5L_lig-to-1C5X_pv.fp

echo 'Ligand formal charge of -2.'
$SCHRODINGER/run $COMBINDHOME/ifp/fp.py -poses 1 -mode pv -input_file inputs/1OWH_lig-to-1C5X_pv.maegz -output_file outputs/1OWH_lig-to-1C5X_pv.fp

echo 'Longest reasonable saltbridge.'
$SCHRODINGER/run $COMBINDHOME/ifp/fp.py -poses 1 -mode pv -input_file inputs/4U16_lig-to-4DAJ_pv.maegz -output_file outputs/4U16_lig-to-4DAJ_pv.fp
