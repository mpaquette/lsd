#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import argparse

from time import time

# import pylab as pl
import nibabel as nib
import numpy as np

from odf_utils import extract_patches




DESCRIPTION = """
Pick the best odf ratio for each voxel based on smoothed AIC in the neighborhood.
"""

EPILOG = """
Michael Paquette, MPI CBS, 2021.
"""


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


def buildArgsParser():

    p = argparse.ArgumentParser(description=DESCRIPTION,
                                epilog=EPILOG,
                                formatter_class=CustomFormatter)

    p.add_argument('--iaic', type=str, nargs='+', default=[],
                             help='Path of the input AICs.')

    p.add_argument('--iodf', type=str, nargs='+', default=[],
                             help='Path of the input ODFs.')

    p.add_argument('--inufo', type=str, nargs='+', default=[],
                             help='Path of the input Nufos.')

    p.add_argument('--ipeak', type=str, nargs='+', default=[],
                             help='Path of the input Peaks.')

    p.add_argument('--ipeaklen', type=str, nargs='+', default=[],
                             help='Path of the input Peak fraction.')

    p.add_argument('--ikMD', type=str, nargs='+', default=[],
                             help='Path of the input kernel MDs.')

    p.add_argument('--mask', type=str, nargs='*', default=[],
                             help='Path of the input mask (one or more).')

    p.add_argument('--ratios', type=float, nargs='+', default=[],
                             help='List of ratios.')  

    p.add_argument('--oodf', type=str,
                             help='Path of the output best ODFs.')

    p.add_argument('--onufo', type=str,
                             help='Path of the output Nufo of the best ODFs.')

    p.add_argument('--opeak', type=str,
                             help='Path of the output Peaks of the best ODFs.')

    p.add_argument('--opeaklen', type=str,
                             help='Path of the output peank length of the best ODFs.')

    p.add_argument('--okMD', type=str,
                             help='Path of the output kernel MD of the best ODFs.')  

    p.add_argument('--oaic', type=str,
                             help='Path of the output best AICs.')

    p.add_argument('--oratio', type=str,
                             help='Path of the output best ratios.')

    return p



def main():
    parser = buildArgsParser()
    args = parser.parse_args()


    # Debug argparse
    # print('Input AIC')
    # for fname in args.iaic:
    #   print(fname)
    # print('Input ODF')
    # for fname in args.iodf:
    #   print(fname)
    # print('Input Mask')
    # for fname in args.mask:
    #   print(fname)
    # print('Input Ratios')
    # for fname in args.ratios:
    #   print(fname)

    # print('Outputs')
    # print(args.oodf)
    # print(args.oaic)
    # print(args.oratio)






    if (args.oodf is None) or (args.oaic is None) or (args.oratio is None):
        print('Need output name(s)')
        return None



    # load and concatenate all the data
    print('Loading AIC data')
    data_img = [nib.load(fname) for fname in args.iaic]
    affine = data_img[0].affine
    data_data = []
    for img in data_img:
        tmp = img.get_fdata()
        # print('data shape = {:}'.format(tmp.shape))
        # need 4D data for the concatenate
        if tmp.ndim == 3:
            tmp = tmp[..., None]
        data_data.append(tmp)
    data_aics = np.concatenate(data_data, axis=3)
    print('Full data shape = {:}'.format(data_aics.shape))
    del data_data


    # load and concatenate all the data
    print('Loading ODF data')
    data_img = [nib.load(fname) for fname in args.iodf]
    # affine = data_img[0].affine
    data_data = []
    for img in data_img:
        tmp = img.get_fdata()
        # print('data shape = {:}'.format(tmp.shape))
        # need 5D data for the concatenate
        if tmp.ndim == 4:
            tmp = tmp[..., None]
        data_data.append(tmp)
    data_odfs = np.concatenate(data_data, axis=4)
    print('Full data shape = {:}'.format(data_odfs.shape))
    del data_data

    # load and concatenate all the data
    print('Loading NUFO data')
    data_img = [nib.load(fname) for fname in args.inufo]
    affine = data_img[0].affine
    data_data = []
    for img in data_img:
        tmp = img.get_fdata()
        # print('data shape = {:}'.format(tmp.shape))
        # need 4D data for the concatenate
        if tmp.ndim == 3:
            tmp = tmp[..., None]
        data_data.append(tmp)
    data_nufos = np.concatenate(data_data, axis=3)
    print('Full data shape = {:}'.format(data_nufos.shape))
    del data_data

    # load and concatenate all the data
    print('Loading PEAK data')
    data_img = [nib.load(fname) for fname in args.ipeak]
    # affine = data_img[0].affine
    data_data = []
    for img in data_img:
        tmp = img.get_fdata()
        # print('data shape = {:}'.format(tmp.shape))
        # need 5D data for the concatenate
        if tmp.ndim == 4:
            tmp = tmp[..., None]
        data_data.append(tmp)
    data_peaks = np.concatenate(data_data, axis=4)
    print('Full data shape = {:}'.format(data_peaks.shape))
    del data_data

    # load and concatenate all the data
    print('Loading PEAKLEN data')
    data_img = [nib.load(fname) for fname in args.ipeaklen]
    # affine = data_img[0].affine
    data_data = []
    for img in data_img:
        tmp = img.get_fdata()
        # print('data shape = {:}'.format(tmp.shape))
        # need 5D data for the concatenate
        if tmp.ndim == 4:
            tmp = tmp[..., None]
        data_data.append(tmp)
    data_peaklens = np.concatenate(data_data, axis=4)
    print('Full data shape = {:}'.format(data_peaklens.shape))
    del data_data

    # load and concatenate all the data
    print('Loading kernel MDs data')
    data_img = [nib.load(fname) for fname in args.ikMD]
    affine = data_img[0].affine
    data_data = []
    for img in data_img:
        tmp = img.get_fdata()
        # print('data shape = {:}'.format(tmp.shape))
        # need 4D data for the concatenate
        if tmp.ndim == 3:
            tmp = tmp[..., None]
        data_data.append(tmp)
    data_kMDs = np.concatenate(data_data, axis=3)
    print('Full data shape = {:}'.format(data_kMDs.shape))
    del data_data


    # load and multiply all the mask
    print('Loading Mask')
    mask = np.ones(data_aics.shape[:3], dtype=np.bool)
    mask_data = [nib.load(fname).get_fdata().astype(np.bool) for fname in args.mask]
    for tmp in mask_data:
        mask = np.logical_and(mask, tmp)
    print('Final mask has {:} voxels ({:.1f} % of total)'.format(mask.sum(), 100*mask.sum()/np.prod(data_aics.shape[:3])))
    del mask_data

    
    ratios = np.array(args.ratios)


    # build 3x3x3 damped kernel
    kn = 3
    padsize = int((kn-1)/2)
    XX,YY,ZZ = np.meshgrid(range(kn), range(kn), range(kn))
    sigma_kernel = 0.5
    kernel = np.exp(-(1/2)*((((XX-(padsize))**2 + (YY-(padsize))**2 + (ZZ-(padsize))**2)/sigma_kernel)**0.5))/ (sigma_kernel*np.sqrt(2*np.pi))
    kernel = kernel[..., None] / kernel.sum() # sum = 1


    # computed kernel aic maps
    dataloop = np.pad(data_aics, ((padsize,), (padsize,), (padsize,), (0,)))    
    maskloop = np.pad(mask, padsize)

    neigh_sum_aic = np.zeros(data_aics.shape)

    data_patches = extract_patches(dataloop, (kn,kn,kn,ratios.shape[0]), (1,1,1,1), flatten=False)
    mask_patches = extract_patches(maskloop, (kn,kn,kn), (1,1,1), flatten=False)

    for xyz in np.ndindex(data_aics.shape[:3]):
        if mask[xyz]:
            datablock = data_patches[xyz+(0,)] * kernel
            maskblock = mask_patches[xyz]

            neigh_sum_aic[xyz] = np.sum(datablock[maskblock], axis=0)





    # find best aic
    print('Find neighboor-kernel best AIC')
    best_idx = np.argmin(neigh_sum_aic, axis=3)


    best_aic = np.zeros(mask.shape)
    best_nufo = np.zeros(mask.shape)
    best_kMD = np.zeros(mask.shape)
    best_ratio = np.zeros(mask.shape)
    best_odf = np.zeros(mask.shape+(data_odfs.shape[3],))
    best_peak = np.zeros(mask.shape+(data_peaks.shape[3],))
    best_peaklen = np.zeros(mask.shape+(data_peaklens.shape[3],))
    # repack everything
    for xyz in np.ndindex(mask.shape):
        best_aic[xyz] = data_aics[xyz][best_idx[xyz]]
        best_nufo[xyz] = data_nufos[xyz][best_idx[xyz]]
        best_kMD[xyz] = data_kMDs[xyz][best_idx[xyz]]
        best_ratio[xyz] = ratios[best_idx[xyz]]
        best_odf[xyz] = data_odfs[xyz][:, best_idx[xyz]]
        best_peak[xyz] = data_peaks[xyz][:, best_idx[xyz]]
        best_peaklen[xyz] = data_peaklens[xyz][:, best_idx[xyz]]



    nib.Nifti1Image(best_aic, affine).to_filename(args.oaic)
    nib.Nifti1Image(best_ratio, affine).to_filename(args.oratio)
    nib.Nifti1Image(best_odf, affine).to_filename(args.oodf)
    nib.Nifti1Image(best_nufo, affine).to_filename(args.onufo)
    nib.Nifti1Image(best_kMD, affine).to_filename(args.okMD)
    nib.Nifti1Image(best_peak, affine).to_filename(args.opeak)
    nib.Nifti1Image(best_peaklen, affine).to_filename(args.opeaklen)



if __name__ == "__main__":
        main()