# Preprocessing

Though some functionality does exist on the graphical interface, this executable is primarily designed to be used from the command line. 

## Functionality offered

- Histogram Matching: needs a reference image and the number of bins and quantiles to match
- Z-Score Normalization
- P1-P2 Normalization: used for DL-based skull stripping
- Bias Correction (N3 and N4)
- Smoothing (Susan)
- Registration (Rigid + Affine + Deformable) using [Greedy](https://sites.google.com/view/greedyreg/about)
- Rescale Images

## To Do

For new features, please put a pull or feature request via github.com/CBICA/CaPTk.
