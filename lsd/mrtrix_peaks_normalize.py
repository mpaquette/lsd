import numpy as np
import nibabel as nib
from tqdm import tqdm
from lsd.odf_utils import vec_angle # in raidans, with proper 180 sym and normalization, and handles 0s



# def main(mrtrix_peaks_fname, config_fname, output_basename):
def main(mrtrix_peaks_fname, peak_rel_th, peak_ang_sep, output_basename):


    PEAK_REL_TH = float(peak_rel_th)
    PEAK_ANG_SEP = float(peak_ang_sep)
    # f = open(config_fname, 'r')
    # config_lines = f.readlines()
    # f.close()
    # lines = [raw_lines for raw_lines in config_lines if raw_lines[0]!='#']
    #
    # dic_idx = [l.strip().split('=')[0] for l in lines]
    # dic_val = [l.strip().split('=')[1] for l in lines]
    # param = dict(zip(dic_idx, dic_val))
    #
    # PEAK_REL_TH = float(param['PEAK_REL_TH']) # fraction
    # PEAK_ANG_SEP = float(param['PEAK_ANG_SEP']) # degrees
    PEAK_ANG_SEP_RAD = PEAK_ANG_SEP*(np.pi/180)


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
    # Compute relative length
    mrtrix_rel_len = mrtrix_len / np.nanmax(mrtrix_len, axis=3)[..., None]
    # Detect small lobe len violations
    mask_too_small_peak = np.logical_and(mrtrix_rel_len < PEAK_REL_TH, ~np.isnan(mrtrix_rel_len))
    # remove too small peaks
    # RECOMPUTE nufo and len
    mrtrix_peaks_reshape[mask_too_small_peak] = np.nan
    # NUFO computed from the number of non NaN entries
    nufo = (~np.isnan(mrtrix_peaks)).sum(axis=3) // 3
    # norm of each peak with nans for non peaks
    mrtrix_len = np.linalg.norm(mrtrix_peaks_reshape, axis=4)

    # detect peaks that are too closed, in order of size
    mask_too_close_peak = np.zeros(mrtrix_peaks_reshape.shape[:4], dtype=bool)
    for idx in tqdm(np.ndindex(mrtrix_peaks_reshape.shape[:3]), total=np.prod(mrtrix_peaks_reshape.shape[:3])):
        voxpeak = mrtrix_peaks_reshape[idx][:nufo[idx]] # Nx3
        for i in range(0, voxpeak.shape[0]-1): # big peak
            if ~mask_too_close_peak[idx + (i,)]:
                for j in range(i+1, voxpeak.shape[0]): # smaller peak
                    if vec_angle(voxpeak[i], voxpeak[j]) < PEAK_ANG_SEP_RAD:
                        mask_too_close_peak[idx + (j,)] = True # kill it

    # remove peaks too close
    # RECOMPUTE nufo and len
    mrtrix_peaks_reshape[mask_too_close_peak] = np.nan
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
    nib.Nifti1Image(nufo.astype(np.int32), affine).to_filename(output_basename + '_nufo.nii.gz')
    nib.Nifti1Image(peak_frac.astype(np.float32), affine).to_filename(output_basename + '_peakfrac.nii.gz')
    nib.Nifti1Image(peak_dir.astype(np.float32), affine).to_filename(output_basename + '_peakdir.nii.gz')


if __name__ == "__main__":
    import sys
    main(*sys.argv[1:])




