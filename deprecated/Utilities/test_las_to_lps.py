#! /usr/bin/env python

import orientations_captk as orn
import numpy as np
from nibabel import load as load_nifti

las_nifti_filename = 'dwi_las.nii.gz'
las_bvec_filename = 'dwi_las.bvec'
las_bval_filename = 'dwi_las.bval' 

lps_bvec_filename = 'dwi_lps.bvec'
lps_bvec_test_filename = 'dwi_lps_TEST.bvec'

# read the dwi
dwi_las_img = load_nifti(las_nifti_filename)
# get the affine (sform/qform) from the image header 
affine = dwi_las_img.affine
# read the bval, bvec 
bvals = np.loadtxt(las_bval_filename)
bvecs = np.loadtxt(las_bvec_filename)
# read the correct lps bvecs for comparison 
bvecs_lps_true = np.loadtxt(lps_bvec_filename)

print('bvals.shape: {}'.format(bvals.shape))
print('bvecs.shape: {}'.format(bvecs.shape))
print('LAS bvecs:')
print(bvecs)

# convert image affine to three-letter tuple 
ornt_code_in = orn.aff2axcodes(affine)
print('Input nifti is {} orientation'.format(''.join(ornt_code_in)))

# desired output orientation is LPS 
ornt_code_out = ('L','P','S')

# do the required reorientation of bvecs from the input and output codes
bvecs_lps = orn.reorient_bvec(bvecs, ornt_code_in, ornt_code_out)
print('LPS bvecs:')
print(bvecs_lps)

# save result 
np.savetxt(lps_bvec_test_filename, bvecs_lps, fmt='%0.8f')

# difference should be zero 
print('Difference between LPS bvecs:')
print(bvecs_lps - bvecs_lps_true)