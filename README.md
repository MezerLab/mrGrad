# mrGrad

This code accompanies two paper that use spatial profiles to map subcortical structures.
The first paper is "Mapping Microstructural Gradients of the Human Striatum in Normal Aging and Parkinson’s Disease" (Drori, Berman, and Mezer, 2021, in press).
This code automatically computes the main axes of the striatum and other subcortical structures in MRI images, based on the ROI voxel coordinates, using SVD on the individual subject level. It then calculates the spatial functions of MRI intensities (or qMRI values) along these axes on the individual and the group levels. Visualization options are also provided for the ROI axes and the MRI functions.
Use mrGrad_run.m script for an example run.
The second paper is "Spatial Profiles of the Midbrain as Sensitive Measures for Underlying Structure" (submitted).
It a script for applying the code to the midbrain is under the directory _MidBrainProfiles_. 
The appropriate example data can be found in each directory.

Requirements:
* MATLAB. Available from: http://www.mathworks.com/products/matlab/
* Vistasoft MATLAB toolbox. Available from: https://github.com/vistalab/vistasoft
Recommended:
* boundedline-pkg MATLAB toolbox. Available from: https://github.com/kakearney/boundedline-pkg
