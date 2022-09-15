#!/usr/bin/env python3

import argparse
import numpy as np
import nibabel as nib
from time import time

from dipy.io import read_bvals_bvecs
# from dipy.core.sphere import HemiSphere
from dipy.core.gradients import gradient_table
# from dipy.sims.voxel import multi_tensor

from odf_utils import true_MD_func, eigenval_from_param, gen_signal_tensor_conv

# from dipy.data import get_sphere
# reg_sphere = get_sphere('symmetric362')

from multiprocessing import cpu_count, Pool

from itertools import combinations as comb


DESCRIPTION = """
compute AIC from all peaks
"""

np.set_printoptions(precision=2)

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
    p.add_argument('--idirs', type=str,
                            help='Name of the input peak dirs')
    p.add_argument('--ilen', type=str,
                            help='Name of the input peak len')
    p.add_argument('--inufo', type=str,
                            help='Name of the input nufo')
    p.add_argument('--sigma', type=str,
                            help='Name of the input sigma map')
    p.add_argument('--ratio', type=float, default = 2., 
                            help='ratio of ODF sharpening')
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
    return np.sum(gaussian_log_likelihood(diffs, sigma))


def main():

    parser = buildArgsParser()
    args = parser.parse_args()

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


    print('Load peak extraction')
    dirs = nib.load(args.idirs).get_fdata()
    lens = nib.load(args.ilen).get_fdata()
    nufo = nib.load(args.inufo).get_fdata()

    # # ghetto number detection
    # if args.sigma.replace('.', '').replace(',','').isdigit():
    #     sigma = float(args.sigma)*np.ones(data.shape[:3])
    # else:
    #     sigma = nib.load(args.sigma).get_fdata()
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
    mask = np.logical_and(mask, sigma > 0)
    # remove nufo = 0 from mask
    mask = np.logical_and(mask, nufo > 0)


    ratio = args.ratio
    NCORE = min(args.cores, cpu_count())



    N_data = data.shape[3]
    N_b0s = gtab.b0s_mask.sum()
    N_dwi = N_data - N_b0s
    N_dirs = dirs.shape[3]//3

    # save the computed kernel MD
    kernel_MDs = np.ones(data.shape[:3])

    # build a function to estimate kernel MD from signal SM
    MD_from_SM = true_MD_func(meanbval=gtab.bvals[~gtab.b0s_mask].mean(), ratio=ratio, minMD=0.01e-3, maxMD=3e-3, N_MD=3000)

    global _aic_loop # this is a hack to make the local function pickleable
    def _aic_loop(data):
        # GTAB is defined outside the function because bad code
        data_vox = data[:N_dwi]
        peak_dir_vox = data[N_dwi:N_dwi+3*N_dirs].reshape((-1, 3))
        peak_len_vox = data[N_dwi+3*N_dirs:N_dwi+4*N_dirs]
        sigma = data[N_dwi+4*N_dirs]
        ratio = data[N_dwi+4*N_dirs+1]

        
        # mean_bval = gtab.bvals[~gtab.b0s_mask].mean()
        # MD_est = np.log(data_vox.mean())/(-mean_bval)
        MD_est = MD_from_SM(data_vox.mean())

        # meval = np.array([1.0, 1/ratio, 1/ratio])
        # K = MD_est / ((1+2/ratio)/3)
        # meval *= K
        Dpar, Dperp = eigenval_from_param(MD_est, ratio)
        # meval = np.array([Dpar, Dperp, Dperp])
        

        vox_npeak = (peak_len_vox>0).sum()
        # use all peaks
        Npeaks = vox_npeak


        # mevals = np.repeat(meval[None, :], Npeaks, axis=0)

        # dirss = peak_dir_vox[:Npeaks]
        # vol_frac = peak_len_vox[:Npeaks]
        # vol_frac /= vol_frac.sum()

        # noiseless_signal, _ = multi_tensor(gtab, mevals=mevals, S0=1, angles=dirss, fractions=100*vol_frac, snr=None)
        noiseless_signal = gen_signal_tensor_conv(gtab, peak_dir_vox[:Npeaks], peak_len_vox[:Npeaks], Dpar, Dperp)



        diffs = data_vox - noiseless_signal[~gtab.b0s_mask]

        loglikelihood = multigaussian_log_likelihood(diffs, sigma)
        criteria_value = aic(loglikelihood, dof=3*Npeaks)

        return criteria_value, MD_est


    print('Concat data')
    concat_data = np.concatenate((data[..., ~gtab.b0s_mask], dirs, lens, sigma[..., None], ratio*np.ones(dirs.shape[:3]+(1,))), axis=3)
    print(concat_data.shape)

    print('Starting AIC all peaks')
    aic_value = np.zeros(concat_data.shape[:3])
    #
    start_time = time()
    loop_output = np.ones(data.shape[:3] + (2,))
    # maybe need chucksize
    with Pool(processes=NCORE) as pool:
        loop_output[mask] = pool.map(_aic_loop, concat_data[mask])
    aic_value = loop_output[..., 0]
    kernel_MDs = loop_output[..., 1]
    end_time = time()
    print('Elapsed time  = {:.2f} s'.format(end_time - start_time))

    # save AIC values
    nib.Nifti1Image(aic_value.astype(np.float), affine).to_filename(args.oaic)
    # save kernel MD
    nib.Nifti1Image(kernel_MDs.astype(np.float), affine).to_filename(args.oMD)


if __name__ == "__main__":
    main()





