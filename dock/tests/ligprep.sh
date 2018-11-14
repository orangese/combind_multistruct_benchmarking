
echo 'pht 2.0 results in incorrect protonation state'
$SCHRODINGER/utilities/prepwizard -WAIT -noepik -noprotassign -noimpref inputs/raw_1F5L.mae prepwizard_1F5L.mae
mv prepwizard_1F5L.mae inputs/

$SCHRODINGER/epik -WAIT -ph 7.0 -pht 2.0 -imae inputs/prepwizard_1F5L.mae -omae inputs/epik_1F5L.mae
$SCHRODINGER/epik -WAIT -ph 7.0 -pht 1.0 -imae inputs/prepwizard_1F5L.mae -omae inputs/epik_1F5L_pht1.mae

mv *.log inputs/
