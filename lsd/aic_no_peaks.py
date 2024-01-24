#!/usr/bin/env python3

import argparse
import numpy as np
import nibabel as nib
from time import time

from dipy.io import read_bvals_bvecs
from dipy.core.gradients import gradient_table

from odf_utils import true_MD_func

# from multiprocessing import cpu_count
from tqdm.contrib.concurrent import process_map

from tqdm import tqdm


DESCRIPTION = """
Fit and tompute AIC for peak-less model
"""

np.set_printoptions(precision=3)

def buildArgsParser():
    p = argparse.ArgumentParser(description=DESCRIPTION)
    p.add_argument('--data', type=str,
                            help='Name of the input noisy dwi')
    p.add_argument('--bval', type=str,
                            help='Name of the input bval')
    p.add_argument('--bvec', type=str,
                            help='Name of the input bvec')
    p.add_argument('--mask', type=str,
                            help='Optional: Name of mask nii file')
    p.add_argument('--sigma', type=str,
                            help='Name of the input sigma map')
    p.add_argument('--oaic', type=str,
                            help='Name of the output aic values')
    p.add_argument('--oMD', type=str,
                            help='Name of the output kernel MD values')
    p.add_argument('--cores', type=int, default = 1, 
                            help='Number of processes')

    return p


def aic(loglikelihood, dof):
    return 2*dof - 2*loglikelihood

def gaussian_log_likelihood(diff, sigma):
    return (-0.5*diff**2/sigma**2) - np.log(sigma*np.sqrt(2*np.pi))

def multigaussian_log_likelihood(diffs, sigma):
    # iid gaussian
    # still works with multiple different sigma, as long as independant
    return np.sum(gaussian_log_likelihood(diffs, sigma))


def main():

    parser = buildArgsParser()
    args = parser.parse_args()

    # cli_string = "--data /data/pt_02101_dMRI/data/007_C_C_NEGRA/preprocessed/20210203_160848_007_C_C_NEGRA_ID11357_1_1_rr_fullproc/diff//odf/lsd_processing/data_norm.nii.gz --bval /data/pt_02101_dMRI/data/007_C_C_NEGRA/preprocessed/20210203_160848_007_C_C_NEGRA_ID11357_1_1_rr_fullproc/diff//odf/lsd_processing/data_norm.bval --bvec /data/pt_02101_dMRI/data/007_C_C_NEGRA/preprocessed/20210203_160848_007_C_C_NEGRA_ID11357_1_1_rr_fullproc/diff//odf/lsd_processing/data_norm.bvec --mask /data/pt_02101_dMRI/data/007_C_C_NEGRA/preprocessed/20210203_160848_007_C_C_NEGRA_ID11357_1_1_rr_fullproc/diff//data_release//mask.nii.gz --sigma /data/pt_02101_dMRI/data/007_C_C_NEGRA/preprocessed/20210203_160848_007_C_C_NEGRA_ID11357_1_1_rr_fullproc/diff//odf/lsd_processing/sigma_norm.nii.gz --oaic /data/pt_02101_dMRI/data/007_C_C_NEGRA/preprocessed/20210203_160848_007_C_C_NEGRA_ID11357_1_1_rr_fullproc/diff//odf/lsd_processing/aic_no_peak.nii.gz --oMD /data/pt_02101_dMRI/data/007_C_C_NEGRA/preprocessed/20210203_160848_007_C_C_NEGRA_ID11357_1_1_rr_fullproc/diff//odf/lsd_processing/MD_no_peak.nii.gz --cores 1"
    # args = parser.parse_args(cli_string.split(" "))

    print('Load data')
    data_img = nib.load(args.data)
    data = data_img.get_fdata()
    affine = data_img.affine

    bvals, bvecs = read_bvals_bvecs(args.bval, args.bvec)
    gtab = gradient_table(bvals, bvecs)

    if args.mask is None:
        mask = np.ones(data.shape[:3], dtype=bool)
    else:
        mask = nib.load(args.mask).get_fdata().astype(bool)

    sigma = nib.load(args.sigma).get_fdata()


    # clean data
    print('Clean data and sigma')
    data = np.clip(data, 0, 1)
    data[np.isnan(data)] = 0
    data[np.isinf(data)] = 0

    # clean sigma
    sigma = np.clip(sigma, 0, np.inf)
    sigma[np.isnan(sigma)] = 0
    sigma[np.isinf(sigma)] = 0

    # remove sigma=0 voxel from mask
    if sigma.ndim == 3:
        mask = np.logical_and(mask, sigma > 0)
    elif sigma.ndim == 4:
        mask = np.logical_and(mask, np.any(sigma > 0, axis=3))

    # NCORE = min(args.cores, cpu_count())

    start_time = time()

    MD_from_SM = true_MD_func(meanbval=gtab.bvals[~gtab.b0s_mask].mean(), ratio=1.0, minMD=0.01e-3, maxMD=3e-3, N_MD=30000)


    MD_est = np.zeros(data.shape[:3])
    for vox in tqdm(np.ndindex(data.shape[:3]), total=np.prod(data.shape[:3])):
        if mask[vox]:
            MD_est[vox] = MD_from_SM(data[vox].mean())


    aic_value = np.zeros(data.shape[:3])
    for vox in tqdm(np.ndindex(data.shape[:3]), total=np.prod(data.shape[:3])):
        if mask[vox]:
            noiseless_signal = np.exp(-gtab.bvals*MD_est[vox])

            diffs = data[vox] - noiseless_signal

            if sigma.ndim == 3:
                sss = sigma[vox]
            elif sigma.ndim == 4:
                sss = sigma[vox][~gtab.b0s_mask]

            loglikelihood = multigaussian_log_likelihood(diffs[~gtab.b0s_mask], sss)
            aic_value[vox] = aic(loglikelihood, dof=0)


    end_time = time()
    print('Elapsed time  = {:.2f} s'.format(end_time - start_time))

    # save AIC values
    nib.Nifti1Image(aic_value.astype(np.float32), affine).to_filename(args.oaic)
    # save MD
    nib.Nifti1Image(MD_est.astype(np.float32), affine).to_filename(args.oMD)


if __name__ == "__main__":
    main()

