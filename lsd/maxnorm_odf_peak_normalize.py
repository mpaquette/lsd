#!/usr/bin/env python3

import argparse
import numpy as np
import nibabel as nib
import glob
import os
from multiprocessing import cpu_count



DESCRIPTION = """
Normalize odf with mrtrix biggest peak lenght
"""

np.set_printoptions(precision=2)

def buildArgsParser():
    p = argparse.ArgumentParser(description=DESCRIPTION)
    p.add_argument('--odf', type=str,
                            help='Filename of the input maxnorm fodf')
    p.add_argument('--mask', type=str,
                            help='Optional: Name of mask nii file')
    p.add_argument('--outPEAK', type=str,
                            help='Name of the output biggest peak.')
    p.add_argument('--outPEAKNORM', type=str,
                            help='Name of the output biggest peak lenght.')
    p.add_argument('--outPEAKNAN', type=str,
                            help='Name of the output map of missing peak.')
    p.add_argument('--outODF', type=str,
                            help='Name of the output fodf sh.')
    p.add_argument('--th', type=float, default=0.25, 
                            help='Peak extraction relative threshold.')
    p.add_argument('--cores', type=int, default=1, 
                            help='Number of processes.')
    return p



def main():

    parser = buildArgsParser()
    args = parser.parse_args()

    odfpath = args.odf
    maskpath = args.mask
    outputpeak = args.outPEAK
    outputpeaknorm = args.outPEAKNORM
    outputnan = args.outPEAKNAN
    outputodfnorm = args.outODF
    rel_th = args.th
    NCORE = min(args.cores, cpu_count())



    # extract biggest peak with mrtrix
    cmd_peakext = \
    'sh2peaks' + ' ' + odfpath \
               + ' ' + outputpeak \
               + ' ' + '-num'       + ' ' + '1'\
               + ' ' + '-threshold' + ' ' + '{:}'.format(rel_th)\
               + ' ' + '-mask'      + ' ' + maskpath\
               + ' ' + '-nthreads'  + ' ' + '{:}'.format(NCORE)\
               + ' ' + '-force'

    os.system(cmd_peakext)



    # load peak
    img_peak = nib.load(outputpeak)
    peakext = img_peak.get_fdata()

    img_mask = nib.load(maskpath)
    mask = img_mask.get_fdata().astype(bool)

    # compute peak len in mask
    peaklen = np.zeros(mask.shape)
    peaklen = np.linalg.norm(peakext, axis=3)

    # save len
    nib.Nifti1Image(peaklen.astype(np.float32), img_mask.affine).to_filename(outputpeaknorm)



    # load odf sh
    img_odf = nib.load(odfpath)
    odf = img_odf.get_fdata()

    # find the few non masked NaNs
    nanmask = np.logical_and(mask, np.isnan(peaklen))
    nib.Nifti1Image(nanmask.astype(np.float32), img_mask.affine).to_filename(outputnan)

    maskpad = np.pad(mask, (1, 1), mode='constant', constant_values=0)
    peaklenpad = np.pad(peaklen, (1, 1), mode='constant', constant_values=0)

    # interpolate len from nonmasked neighboor
    X,Y,Z = np.where(nanmask)
    print('Found {:} non peak ext vox'.format(X.shape[0]))
    for i in range(X.shape[0]):
        # tmp_mask = mask[X[i]-1:X[i]+2, Y[i]-1:Y[i]+2, Z[i]-1:Z[i]+2]
        tmp_mask = maskpad[X[i]-1+1:X[i]+2+1, Y[i]-1+1:Y[i]+2+1, Z[i]-1+1:Z[i]+2+1]
        tmp_mask[1,1,1] = 0 # remove vox itself
        count = tmp_mask.sum()

        if count > 0:
            # peaklen[X[i], Y[i], Z[i]] = peaklen[X[i]-1:X[i]+2, Y[i]-1:Y[i]+2, Z[i]-1:Z[i]+2][tmp_mask].mean()
            peaklen[X[i], Y[i], Z[i]] = peaklenpad[X[i]-1+1:X[i]+2+1, Y[i]-1+1:Y[i]+2+1, Z[i]-1+1:Z[i]+2+1][tmp_mask].mean()
        else:
            # if theres no non masked neighbor, keep odf as is
            peaklen[X[i], Y[i], Z[i]] = 1.0

        print('Vox ({:}, {:}, {:}) new len is {:.3f}'.format(X[i], Y[i], Z[i], peaklen[X[i], Y[i], Z[i]]))


    # normalize with fixed peak len
    odf = odf / peaklen[..., None]
    odf[~mask] = 0 # null the outside

    # save normalized odf
    print('Saving')
    nib.Nifti1Image(odf.astype(np.float32), img_mask.affine).to_filename(outputodfnorm)



if __name__ == "__main__":
    main()
