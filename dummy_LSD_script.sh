CONFIGFILE=$1

. $CONFIGFILE


starttime=`date +%s`

# TODO:
# basic input "existence" sanity checks

echo $DWI_PATH;
echo $BVAL_PATH;
echo $BVEC_PATH;
echo $MASK_PATH;
echo $SIGMA_PATH;
echo $N_PATH;
echo $OUTPUT_BASE;
echo $RATIOS;
echo $SHMAX;
echo $PEAK_REL_TH;
echo $PEAK_ANG_SEP;  
echo $PATCHSIZE;
echo $NCORE;


ALLCORE=$(nproc --all)
if [[ $NCORE -gt $ALLCORE ]]
then
    NCORE=$ALLCORE;
fi

if [[ $NCORE -eq -1 ]]
then
    NCORE=$ALLCORE;
fi


OUTPUT_FOLDER=$(dirname $OUTPUT_BASE);
PROC_FOLDER=$OUTPUT_FOLDER'/lsd_processing/'
mkdir $PROC_FOLDER



echo 'Make sure Noisemap are 3D'
python3 lsd/make_noisemap.py \
    $SIGMA_PATH \
    $N_PATH \
    $MASK_PATH \
    $PROC_FOLDER'sigma_gen.nii.gz' \
    $PROC_FOLDER'Ns_gen.nii.gz' 

# overwrite paths
SIGMA_PATH=$PROC_FOLDER'sigma_gen.nii.gz'
N_PATH=$PROC_FOLDER'Ns_gen.nii.gz'


####################################
# Normalize Data with b0

echo 'Normalize Data with b0'
python3 lsd/normalize_data.py \
    --in $DWI_PATH \
    --in_sigma $SIGMA_PATH \
    --in_N $N_PATH \
    --mask $MASK_PATH \
    --bval $BVAL_PATH \
    --bvec $BVEC_PATH \
    --out_folder $PROC_FOLDER
#
##################


echo 'Fit CSA odf'
python3 lsd/fit_csa.py \
    $PROC_FOLDER'data_norm.nii.gz' \
    $PROC_FOLDER'data_norm.bval' \
    $PROC_FOLDER'data_norm.bvec' \
    $MASK_PATH \
    $PROC_FOLDER'csa.nii.gz' \
    $NCORE \
    1e-5 \
    0.006 \
    $SHMAX


for RATIO in ${RATIOS[@]};
do
    echo 'Ratio '$RATIO

    if [[ ! -e $PROC_FOLDER'kernelMD_csa_sharp_r'$RATIO'.nii.gz' ]]
    then

        # Descoteaux "sharpening odf" deconv
        python3 lsd/sharpen_sh_parallel.py \
                --in $PROC_FOLDER'csa.nii.gz' \
                --out $PROC_FOLDER'csa_sharp_r'$RATIO'.nii.gz' \
                --mask $MASK_PATH \
                --ratio $RATIO \
                --tau 0.1 \
                --lambda 1. \
                --csa_norm True \
                --cores $NCORE

        # Normalize ODF with max=1
        # This allows for absolute thresholds to behave like relative thresholds
        python3 lsd/sh_odf_normalize.py \
                $PROC_FOLDER'csa_sharp_r'$RATIO'.nii.gz' \
                $PROC_FOLDER'csa_sharp_r'$RATIO'_norm.nii.gz'

        # Non linear peak extraction
        sh2peaks $PROC_FOLDER'csa_sharp_r'$RATIO'_norm.nii.gz' \
                 $PROC_FOLDER'csa_sharp_r'$RATIO'_norm_peakext.nii.gz' \
                 -num 10 \
                 -threshold $PEAK_REL_TH \
                 -mask $MASK_PATH \
                 -nthreads $NCORE \
                 -force

        # Compute normalize directions and peak fractions
        python3 lsd/mrtrix_peaks_normalize.py \
                $PROC_FOLDER'csa_sharp_r'$RATIO'_norm_peakext.nii.gz' \
                $PEAK_REL_TH \
                $PEAK_ANG_SEP \
                $PROC_FOLDER'csa_sharp_r'$RATIO

        # Compute AIC
        python3 lsd/compute_aic_all_peaks.py \
                --data $PROC_FOLDER'data_norm.nii.gz' \
                --bval $PROC_FOLDER'data_norm.bval' \
                --bvec $PROC_FOLDER'data_norm.bvec' \
                --mask $MASK_PATH \
                --inufo $PROC_FOLDER'csa_sharp_r'$RATIO'_nufo.nii.gz' \
                --idirs $PROC_FOLDER'csa_sharp_r'$RATIO'_peakdir.nii.gz' \
                --ilen $PROC_FOLDER'csa_sharp_r'$RATIO'_peakfrac.nii.gz' \
                --sigma $PROC_FOLDER'sigma_norm.nii.gz' \
                --ratio $RATIO \
                --oaic $PROC_FOLDER'aic_csa_sharp_r'$RATIO'.nii.gz' \
                --oMD $PROC_FOLDER'kernelMD_csa_sharp_r'$RATIO'.nii.gz' \
                --cores $NCORE
    else
        echo 'Output files exist, skipping';
    fi
done



# stack odf and aic filename in a list, in order of increasing ratios
declare -a ODFFILELIST
declare -a AICFILELIST
declare -a NUFOFILELIST
declare -a PEAKFILELIST
declare -a PEAKLENFILELIST
declare -a KERNELMDFILELIST
for RATIO in ${RATIOS[@]};
do
    TMP_ODFFILE=$PROC_FOLDER'csa_sharp_r'$RATIO'.nii.gz';
    ODFFILELIST[${#ODFFILELIST[@]}+1]=$TMP_ODFFILE;
    #
    TMP_AICFILE=$PROC_FOLDER'aic_csa_sharp_r'$RATIO'.nii.gz';
    AICFILELIST[${#AICFILELIST[@]}+1]=$TMP_AICFILE;
    #
    TMP_NUFOFILE=$PROC_FOLDER'csa_sharp_r'$RATIO'_nufo.nii.gz';
    NUFOFILELIST[${#NUFOFILELIST[@]}+1]=$TMP_NUFOFILE;
    #
    TMP_PEAKFILE=$PROC_FOLDER'csa_sharp_r'$RATIO'_peakdir.nii.gz';
    PEAKFILELIST[${#PEAKFILELIST[@]}+1]=$TMP_PEAKFILE;
    #
    TMP_PEAKLENFILE=$PROC_FOLDER'csa_sharp_r'$RATIO'_peakfrac.nii.gz';
    PEAKLENFILELIST[${#PEAKLENFILELIST[@]}+1]=$TMP_PEAKLENFILE;
    #
    TMP_KERNELMDFILE=$PROC_FOLDER'kernelMD_csa_sharp_r'$RATIO'.nii.gz';
    KERNELMDFILELIST[${#KERNELMDFILELIST[@]}+1]=$TMP_KERNELMDFILE;
done


if [[ $PATCHSIZE -gt 1 ]]
then
    # Picks the ratio with lowest AIC for each voxel in neighborhood
    python3 lsd/combine_aic_neigh.py \
            --iaic ${AICFILELIST[@]} \
            --iodf ${ODFFILELIST[@]} \
            --inufo ${NUFOFILELIST[@]} \
            --ipeak ${PEAKFILELIST[@]} \
            --ipeaklen ${PEAKLENFILELIST[@]} \
            --ikMD ${KERNELMDFILELIST[@]} \
            --mask $MASK_PATH \
            --ratios ${RATIOS[@]} \
            --oodf $OUTPUT_FOLDER'/odf_best_aic.nii.gz' \
            --onufo $OUTPUT_FOLDER'/nufo_best_aic.nii.gz' \
            --opeak $OUTPUT_FOLDER'/peak_best_aic.nii.gz' \
            --opeaklen $OUTPUT_FOLDER'/peaklen_best_aic.nii.gz' \
            --okMD $OUTPUT_FOLDER'/kernel_MD_best_aic.nii.gz' \
            --oaic $OUTPUT_FOLDER'/aic_best_aic.nii.gz' \
            --oratio $OUTPUT_FOLDER'/ratio_best_aic.nii.gz';
else
    # Picks the ratio with lowest AIC for each voxel
    python3 lsd/combine_aic.py \
            --iaic ${AICFILELIST[@]} \
            --iodf ${ODFFILELIST[@]} \
            --inufo ${NUFOFILELIST[@]} \
            --ipeak ${PEAKFILELIST[@]} \
            --ipeaklen ${PEAKLENFILELIST[@]} \
            --ikMD ${KERNELMDFILELIST[@]} \
            --mask $MASK_PATH \
            --ratios ${RATIOS[@]} \
            --oodf $OUTPUT_FOLDER'/odf_best_aic.nii.gz' \
            --onufo $OUTPUT_FOLDER'/nufo_best_aic.nii.gz' \
            --opeak $OUTPUT_FOLDER'/peak_best_aic.nii.gz' \
            --opeaklen $OUTPUT_FOLDER'/peaklen_best_aic.nii.gz' \
            --okMD $OUTPUT_FOLDER'/kernel_MD_best_aic.nii.gz' \
            --oaic $OUTPUT_FOLDER'/aic_best_aic.nii.gz' \
            --oratio $OUTPUT_FOLDER'/ratio_best_aic.nii.gz';
fi


# Normalize ODF with max=1
# This allows for absolute thresholds to behave like relative thresholds
python3 lsd/sh_odf_normalize.py \
        $OUTPUT_FOLDER'/odf_best_aic.nii.gz' \
        $OUTPUT_FOLDER'/odf_best_aic_maxnorm.nii.gz'


endtime=`date +%s`
runtime=$((endtime-starttime))

echo 'execution time = '$runtime' seconds'


