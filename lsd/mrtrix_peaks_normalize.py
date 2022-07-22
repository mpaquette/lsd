import numpy as np
import nibabel as nib
# from time import time


def main(mrtrix_peaks_fname, output_basename):

    # Load Mrtrix sh2peaks output
    # it contains the peak in format (X, Y, Z, 3*NPEAKS)
    # with nans as fillers
    # peak lenghts are equivalent (but not equal) to odf lobes
    img = nib.load(mrtrix_peaks_fname)
    mrtrix_peaks = img.get_fdata()
    affine = img.affine

    # reshape into (X, Y, Z, NPEAKS, 3)
    mrtrix_peaks_reshape = mrtrix_peaks.reshape(mrtrix_peaks.shape[:3]+(mrtrix_peaks.shape[3]//3, 3))

    # NUFO computed from the number of non NaN entries
    nufo = (~np.isnan(mrtrix_peaks)).sum(axis=3) // 3

    # norm of each peak with nans for non peaks
    mrtrix_len = np.linalg.norm(mrtrix_peaks_reshape, axis=4)

    # compute fractions and replace nans by 0
    peak_frac = mrtrix_len / np.nansum(mrtrix_len, axis=3)[..., None]
    peak_frac[np.isnan(peak_frac)] = 0

    # dir has norm 1 vectors for peaks and zeros elsewhere
    peak_dir = mrtrix_peaks_reshape / mrtrix_len[..., None]
    peak_dir[np.isnan(peak_dir)] = 0
    peak_dir = peak_dir.reshape(peak_dir.shape[:3] + (peak_dir.shape[3]*3, ))

    # save
    nib.Nifti1Image(nufo.astype(np.int), affine).to_filename(output_basename + '_nufo.nii.gz')
    nib.Nifti1Image(peak_frac.astype(np.float32), affine).to_filename(output_basename + '_peakfrac.nii.gz')
    nib.Nifti1Image(peak_dir.astype(np.float32), affine).to_filename(output_basename + '_peakdir.nii.gz')


if __name__ == "__main__":
    import sys
    main(*sys.argv[1:])




