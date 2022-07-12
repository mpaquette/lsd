#!/usr/bin/env python3

import argparse
import numpy as np
import nibabel as nib


DESCRIPTION = """
Convert Sigma and Number of Coil maps into proper format
"""

def buildArgsParser():
    p = argparse.ArgumentParser(description=DESCRIPTION)

    p.add_argument('sigma', type=str,
                            help='Path to sigma map or single sigma value')
    p.add_argument('Ns', type=str,
                            help='Path to number of coil map or single N value')
    p.add_argument('ref', type=str,
                            help='Reference image')
    p.add_argument('osigma', type=str,
                            help='Path to sigma map output')
    p.add_argument('oNs', type=str,
                            help='Path to number of coil map output')

    return p

def main():

    parser = buildArgsParser()
    args = parser.parse_args()

    ref_img = nib.load(args.ref)

    if args.sigma.replace('.', '').replace(',','').isdigit():
        sigma = float(args.sigma)*np.ones(ref_img.shape[:3])
    else:
        sigma = nib.load(args.sigma).get_fdata()

    if args.Ns.replace('.', '').replace(',','').isdigit():
        Ns = float(args.Ns)*np.ones(ref_img.shape[:3])
    else:
        Ns = nib.load(args.Ns).get_fdata()


    nib.Nifti1Image(sigma.astype(np.float), ref_img.affine).to_filename(args.osigma)
    nib.Nifti1Image(Ns.astype(np.float), ref_img.affine).to_filename(args.oNs)



if __name__ == "__main__":
    main()

