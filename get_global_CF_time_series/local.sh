echo "You need to have an activate conad env that has access to all needed packages"
echo "Consider 'conda activate geo_stuff' based on README file: ../get_regional_CF_time_series_on_linux-mac/README.md"

for YEAR in {1980..2019}; do
    echo ""
    echo ""
    echo ""
    echo ""
    echo ""
    echo ""
    echo ""
    date
    echo ""
    echo $YEAR
    source step0_get_temperature_profile.sh 
    date
done
