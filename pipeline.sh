COMBIND_CODE="/scratch/PI/rondror/augustine/combind_testing/augustines-combind"
WRITE_DIR="/scratch/PI/rondror/combind/nonbinders"
PROTEINS='B1AR'
PREPARATION_STAGE='2'
# Stages:
# 2: Get ligands and process ligands
# c: Redo chembl/ files
# 3: Pick helpers
# 4: Dock helpers
# 5: Compute fingerprints
# 6: Compute mcss

if [[ $COMBIND_HAS_BEEN_SETUP != "OUI" ]]
then
    COMBIND_HAS_BEEN_SETUP="OUI"
    cd $COMBIND_CODE
    source $COMBIND_CODE/setup_combind.sh
    cd $COMBIND_CODE
fi
cd $WRITE_DIR



# To prepare ligands
$SCHRODINGER/run $COMBIND_CODE/main.py prepare $PREPARATION_STAGE $PROTEINS

# To compute dude statistics
# $SCHRODINGER/run /scratch/PI/rondror/augustine/combind_testing/augustines-combind/main.py dude_statistics $PROTEINS

cd $COMBIND_CODE
