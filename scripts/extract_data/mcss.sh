
prot=$1
for lig1 in $(ls $PI_SCRATCH/combind/bpp_data/$prot/ligands/prepared_ligands/); do
    if [ ${lig1:0:6} != 'CHEMBL' ]; then
	for lig2 in $(ls $PI_SCRATCH/combind/bpp_data/$prot/ligands/prepared_ligands/); do
	    if [ ${lig2:0:6} != 'CHEMBL' ] && [ $lig1 \< $lig2 ]; then
		echo $lig1 $lig2
		$SCHRODINGER/run visualize_mcss.py $lig1 $lig2 $prot
	    fi
	done
    fi
done
