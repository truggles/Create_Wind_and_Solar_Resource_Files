

DATE=20210113v5


# Apparently I need to deactivate my conda env to allow
# SLURM to pick up the new env on the other side...
conda deactivate

echo "If having trouble with that conad envs, make sure that you deactivate the base environment too!"

# NYIS, ERCO, PJM 2018 files can use `selected_mask_outfile_US_BAS.nc`
# FR can use `selected_mask_outfile_FR.nc`
# from here: https://github.com/truggles/Inter-annual_Variability_Residual_Load/tree/main/shape_files

        #"TI" \
        #"FR" \
        #"NYIS_2018" \
        #"ERCO_2018" \
        #"PJM_2018" \

for REGION in \
        "TEX" \
        "NYS" \
        ; do

    for METHOD in {1..3}; do
    #for METHOD in {1..1}; do

        for YEAR in {1980..2019}; do

            for CFS in "scf" "wcf"; do
            #for CFS in "temp"; do
                echo "${DATE} ${REGION} METHOD:${METHOD} YEAR:${YEAR} cfs:${CFS}"
                export EXTRA_ARGS="${YEAR} ${REGION} ${METHOD} ${DATE} ${CFS}"
                export SBATCH_JOB_NAME=cfs_${DATE}_${REGION}_mthd${METHOD}_y${YEAR}_${CFS}
                sbatch step2_run.sh
            done

        done

    done

done
