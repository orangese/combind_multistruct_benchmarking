# Downloaded PDB: 1L2S. Split into protein, ligand in Maestro.

# Download known binders from ChEMBL 27
python ~/combind/scripts/chembl.py CHEMBL2026 --homologous --activity_type IC50
python ~/combind/scripts/split_chembl.py CHEMBL2026_IC50.csv

combind ligprep CHEMBL2026_IC50_nM.smi ligands
combind ligprep ampc_screen.smi ligands --screen

combind dock structures/grids/1L2S/1L2S.zip docking ligands/CHEMBL122450/CHEMBL122450.maegz ligands/CHEMBL124416/CHEMBL124416.maegz ligands/CHEMBL127782/CHEMBL127782.maegz ligands/CHEMBL263746/CHEMBL263746.maegz ligands/CHEMBL295322/CHEMBL295322.maegz ligands/CHEMBL297805/CHEMBL297805.maegz ligands/CHEMBL331090/CHEMBL331090.maegz ligands/CHEMBL340707/CHEMBL340707.maegz ligands/CHEMBL416447/CHEMBL416447.maegz ligands/CHEMBL42528/CHEMBL42528.maegz ligands/CHEMBL44813/CHEMBL44813.maegz ligands/CHEMBL44932/CHEMBL44932.maegz ligands/CHEMBL47198/CHEMBL47198.maegz
combind dock structures/grids/1L2S/1L2S.zip docking ligands/ampc_screen/ampc_screen.maegz --screen


combind ligprep dopamine_iuphar.smi ligands
combind ligprep dopamine_screen.smi ligands --screen

combind dock structures/grids/5WIU/5WIU.zip docking ligands/dopamine_screen/dopamine_screen.maegz --screen
combind dock structures/grids/5WIU/5WIU.zip docking ligands/A-381393/A-381393.maegz ligands/aripiprazole/aripiprazole.maegz ligands/benperidol/benperidol.maegz ligands/chlorpromazine/chlorpromazine.maegz ligands/clozapine/clozapine.maegz ligands/CP-293019/CP-293019.maegz ligands/eticlopride/eticlopride.maegz ligands/FAUC213/FAUC213.maegz ligands/haloperidol/haloperidol.maegz ligands/L741742/L741742.maegz ligands/L745870/L745870.maegz ligands/L-750667/L-750667.maegz ligands/loxapine/loxapine.maegz ligands/minus-sulpiride/minus-sulpiride.maegz ligands/ML398/ML398.maegz ligands/nemonapride/nemonapride.maegz ligands/NGD_94-1/NGD_94-1.maegz ligands/N-methylspiperone/N-methylspiperone.maegz ligands/perospirone/perospirone.maegz ligands/RBI257/RBI257.maegz ligands/Ro_10-4548/Ro_10-4548.maegz ligands/sertindole/sertindole.maegz ligands/sonepiprazole/sonepiprazole.maegz ligands/spiperone/spiperone.maegz ligands/terguride/terguride.maegz ligands/trifluoperazine/trifluoperazine.maegz ligands/U101958/U101958.maegz ligands/zotepine/zotepine.maegz

combind featurize . docking/A-381393-to-5WIU/A-381393-to-5WIU_pv.maegz docking/aripiprazole-to-5WIU/aripiprazole-to-5WIU_pv.maegz docking/benperidol-to-5WIU/benperidol-to-5WIU_pv.maegz docking/chlorpromazine-to-5WIU/chlorpromazine-to-5WIU_pv.maegz docking/clozapine-to-5WIU/clozapine-to-5WIU_pv.maegz docking/CP-293019-to-5WIU/CP-293019-to-5WIU_pv.maegz docking/eticlopride-to-5WIU/eticlopride-to-5WIU_pv.maegz docking/FAUC213-to-5WIU/FAUC213-to-5WIU_pv.maegz docking/haloperidol-to-5WIU/haloperidol-to-5WIU_pv.maegz docking/L741742-to-5WIU/L741742-to-5WIU_pv.maegz docking/L745870-to-5WIU/L745870-to-5WIU_native_pv.maegz docking/L-750667-to-5WIU/L-750667-to-5WIU_pv.maegz docking/loxapine-to-5WIU/loxapine-to-5WIU_pv.maegz docking/minus-sulpiride-to-5WIU/minus-sulpiride-to-5WIU_pv.maegz docking/ML398-to-5WIU/ML398-to-5WIU_pv.maegz docking/nemonapride-to-5WIU/nemonapride-to-5WIU_native_pv.maegz docking/NGD_94-1-to-5WIU/NGD_94-1-to-5WIU_pv.maegz docking/N-methylspiperone-to-5WIU/N-methylspiperone-to-5WIU_pv.maegz docking/perospirone-to-5WIU/perospirone-to-5WIU_pv.maegz docking/RBI257-to-5WIU/RBI257-to-5WIU_pv.maegz docking/Ro_10-4548-to-5WIU/Ro_10-4548-to-5WIU_pv.maegz docking/sertindole-to-5WIU/sertindole-to-5WIU_pv.maegz docking/sonepiprazole-to-5WIU/sonepiprazole-to-5WIU_pv.maegz docking/spiperone-to-5WIU/spiperone-to-5WIU_pv.maegz docking/terguride-to-5WIU/terguride-to-5WIU_pv.maegz docking/trifluoperazine-to-5WIU/trifluoperazine-to-5WIU_pv.maegz docking/U101958-to-5WIU/U101958-to-5WIU_pv.maegz docking/zotepine-to-5WIU/zotepine-to-5WIU_pv.maegz

combind pose-prediction bpp dopamine_iuphar_native.csv A-381393-to-5WIU_pv aripiprazole-to-5WIU_pv benperidol-to-5WIU_pv chlorpromazine-to-5WIU_pv clozapine-to-5WIU_pv CP-293019-to-5WIU_pv eticlopride-to-5WIU_pv FAUC213-to-5WIU_pv haloperidol-to-5WIU_pv L741742-to-5WIU_pv L745870-to-5WIU_native_pv L-750667-to-5WIU_pv loxapine-to-5WIU_pv minus-sulpiride-to-5WIU_pv ML398-to-5WIU_pv nemonapride-to-5WIU_native_pv NGD_94-1-to-5WIU_pv N-methylspiperone-to-5WIU_pv perospirone-to-5WIU_pv RBI257-to-5WIU_pv Ro_10-4548-to-5WIU_pv sertindole-to-5WIU_pv sonepiprazole-to-5WIU_pv spiperone-to-5WIU_pv terguride-to-5WIU_pv trifluoperazine-to-5WIU_pv U101958-to-5WIU_pv zotepine-to-5WIU_pv
combind extract-top-poses dopamine_iuphar_native.csv bpp/docking



combind featurize screen docking/dopamine_screen-to-5WIU/dopamine_screen-to-5WIU_pv.maegz ../dopamine_iuphar_pv.maegz --no-mcss --screen
combind screen screen/screen.npy screen/screen/gscore/dopamine_screen-to-5WIU_pv.npy --ifp-fname screen/screen/ifp-pair/{}-dopamine_screen-to-5WIU_pv-and-dopamine_iuphar_pv.npy --features hbond,contact,saltbridge
combind apply-scores screen/docking/dopamine_screen-to-5WIU/dopamine_screen-to-5WIU_pv.maegz screen/screen.npy screen/screen_pv.maegz
$SCHRODINGER/utilities/glide_sort -best_by_title -use_prop_d r_i_combind_score  -o screen/screen_combind_pv.maegz screen/screen_pv.maegz
$SCHRODINGER/utilities/glide_sort -best_by_title  -o screen/screen_glide_pv.maegz screen/screen_pv.maegz

combind scores-to-csv screen/screen_glide_pv.maegz screen/glide.csv
combind scores-to-csv screen/screen_combind_pv.maegz screen/combind.csv


combind ligprep melatonin_iuphar.smi ligands
combind ligprep melatonin_screen.smi ligands --screen

combind dock structures/grids/6ME2/6ME2.zip docking ligands/melatonin_screen/melatonin_screen.maegz --screen
combind dock structures/grids/6ME2/6ME2.zip docking ligands/2-iodo-melatonin/2-iodo-melatonin.maegz ligands/2-methoxy-alpha-beta-didehydro-agomelatine/2-methoxy-alpha-beta-didehydro-agomelatine.maegz ligands/5-HEAT/5-HEAT.maegz ligands/6-Cl-MLT/6-Cl-MLT.maegz ligands/6-hydroxymelatonin/6-hydroxymelatonin.maegz ligands/AAE-M-PBP-amine/AAE-M-PBP-amine.maegz ligands/agomelatine/agomelatine.maegz ligands/CBOBNEA/CBOBNEA.maegz ligands/difluoroagomelatine/difluoroagomelatine.maegz ligands/EFPPEA/EFPPEA.maegz ligands/GR_128107/GR_128107.maegz ligands/GR_196429/GR_196429.maegz ligands/IIK7/IIK7.maegz ligands/ISD6/ISD6.maegz ligands/LY_156735/LY_156735.maegz ligands/melatonin/melatonin.maegz ligands/ramelteon/ramelteon.maegz ligands/S24014/S24014.maegz ligands/S24773/S24773.maegz ligands/S26284/S26284.maegz ligands/tasimelteon/tasimelteon.maegz ligands/UCM_793/UCM_793.maegz

combind featurize . docking/2-iodo-melatonin-to-6ME2/2-iodo-melatonin-to-6ME2_native_pv.maegz docking/2-methoxy-alpha-beta-didehydro-agomelatine-to-6ME2/2-methoxy-alpha-beta-didehydro-agomelatine-to-6ME2_pv.maegz docking/5-HEAT-to-6ME2/5-HEAT-to-6ME2_pv.maegz docking/6-Cl-MLT-to-6ME2/6-Cl-MLT-to-6ME2_pv.maegz docking/6-hydroxymelatonin-to-6ME2/6-hydroxymelatonin-to-6ME2_pv.maegz docking/AAE-M-PBP-amine-to-6ME2/AAE-M-PBP-amine-to-6ME2_pv.maegz docking/agomelatine-to-6ME2/agomelatine-to-6ME2_native_pv.maegz docking/CBOBNEA-to-6ME2/CBOBNEA-to-6ME2_pv.maegz docking/difluoroagomelatine-to-6ME2/difluoroagomelatine-to-6ME2_pv.maegz docking/EFPPEA-to-6ME2/EFPPEA-to-6ME2_pv.maegz docking/GR_128107-to-6ME2/GR_128107-to-6ME2_pv.maegz docking/GR_196429-to-6ME2/GR_196429-to-6ME2_pv.maegz docking/IIK7-to-6ME2/IIK7-to-6ME2_pv.maegz docking/ISD6-to-6ME2/ISD6-to-6ME2_pv.maegz docking/LY_156735-to-6ME2/LY_156735-to-6ME2_pv.maegz docking/melatonin-to-6ME2/melatonin-to-6ME2_pv.maegz docking/ramelteon-to-6ME2/ramelteon-to-6ME2_native_pv.maegz docking/S24014-to-6ME2/S24014-to-6ME2_pv.maegz docking/S24773-to-6ME2/S24773-to-6ME2_pv.maegz docking/S26284-to-6ME2/S26284-to-6ME2_pv.maegz docking/tasimelteon-to-6ME2/tasimelteon-to-6ME2_pv.maegz docking/UCM_793-to-6ME2/UCM_793-to-6ME2_pv.maegz

combind pose-prediction . binders.csv 2-iodo-melatonin-to-6ME2_native_pv 6-hydroxymelatonin-to-6ME2_pv difluoroagomelatine-to-6ME2_pv IIK7-to-6ME2_pv melatonin-to-6ME2_pv S26284-to-6ME2_pv 2-methoxy-alpha-beta-didehydro-agomelatine-to-6ME2_pv AAE-M-PBP-amine-to-6ME2_pv EFPPEA-to-6ME2_pv ISD6-to-6ME2_pv ramelteon-to-6ME2_native_pv tasimelteon-to-6ME2_pv 5-HEAT-to-6ME2_pv agomelatine-to-6ME2_native_pv GR_128107-to-6ME2_pv LY_156735-to-6ME2_pv S24014-to-6ME2_pv UCM_793-to-6ME2_pv 6-Cl-MLT-to-6ME2_pv CBOBNEA-to-6ME2_pv GR_196429-to-6ME2_pv S24773-to-6ME2_pv --xtal 2-iodo-melatonin-to-6ME2_native_pv --xtal ramelteon-to-6ME2_native_pv --xtal agomelatine-to-6ME2_native_pv --alpha 1.0 --gc50 -8.0 --features shape,mcss,hbond,saltbridge,contact

combind extract-top-poses binders.csv docking
combind featurize screen screen/docking/melatonin_screen-to-6ME2/melatonin_screen-to-6ME2_pv.maegz bpp/binders_pv.maegz --no-mcss --screen

combind screen screen/screen.npy screen/gscore/melatonin_screen-to-6ME2_pv.npy --ifp-fname screen/ifp-pair/{}-melatonin_screen-to-6ME2_pv-and-binders_pv.npy --shape-fname screen/shape/shape-melatonin_screen-to-6ME2_pv-and-binders_pv.npy
combind apply-scores ../../../subset/docking/subset-to-XTAL/subset-to-XTAL_pv.maegz screen/screen.npy screen/screen_pv.maegz

$SCHRODINGER/utilities/glide_sort -best_by_title -use_prop_d r_i_combind_score  -o screen/screen_combind_pv.maegz screen/screen_pv.maegz
$SCHRODINGER/utilities/glide_sort -best_by_title  -o screen/screen_glide_pv.maegz screen/screen_pv.maegz

combind scores-to-csv screen/screen_glide_pv.maegz screen/glide.csv
combind scores-to-csv screen/screen_combind_pv.maegz screen/combind.csv
