# mrGrad: Deep Brain MRI Gradient Analysis Toolbox

This toolbox accompanies two research papers that utilize spatial profiles (gradients) to map subcortical brain structures using MRI data.

## Paper 1  
**Mapping Microstructural Gradients of the Human Striatum in Normal Aging and Parkinson’s Disease**
*Drori, Berman, and Mezer. Science Advances, 2022*

This code automatically computes the principal axes of the striatum (caudate and putamen) and possibly other subcortical structures (e.g., substantia nigra) in MRI images. It uses **singular value decomposition (SVD)** on ROI voxel coordinates at the individual subject level.

It then calculates **spatial functions** of MRI intensities (or qMRI values) along these axes, both at the individual and group levels.

It outputs **summary results** of the image value spatial profiles and the axes information, as well as **subject-specific NIfTI files of the resulting axis-based segmentations**

Visualization tools are provided for:
- ROI axes
- MRI-derived spatial functions

**Example usage**: Run the script `mrGrad_run.m`

## Paper 2  
**Spatial profiles provide sensitive MRI measures of the midbrain micro- and macrostructure**
*Berman, Drori, and Mezer. NeuroImage, 2022*

A separate script for applying this method to the midbrain is available in the `MidBrainProfiles` directory.
**Example data** for each analysis is included in their respective directories.

## Requirements:
`mrGrad` is a fully self-contained **MATLAB-based** toolbox.

**Required:**
- [MATLAB](http://www.mathworks.com/products/matlab/)


**Recommended:**
- [`boundedline-pkg`](https://github.com/kakearney/boundedline-pkg): A MATLAB toolbox for plotting bounded lines (useful for visualization)
