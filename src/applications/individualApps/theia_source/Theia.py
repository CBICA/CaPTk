import argparse
import sys
import os

#from visualizer 

def verify_type(file):
    ext = os.path.basename(file).split(os.extsep, 1)
    if ext[1] != 'nii.gz':
        parser.error("File doesn't end with 'nii.gz'. Found: {}".format(ext[1]))
    return file


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Reads Nii.gz Files and renders them in 3D.')
    parser.add_argument('-i', type=lambda fn: verify_type(fn), help='an mri scan (nii.gz)')
    parser.add_argument('-m', type=lambda fn: verify_type(fn), help='the segmentation mask (nii.gz)')
    args = parser.parse_args()
    
    os.system(os.getcwd() + "/visualizer/brain_tumor_3d.py -i " + args.i + " -m " + args.m)