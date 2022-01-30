# mrGrad

This code accompnies the paper "Mapping Microstructural Gradients of the Human Striatum in Normal Aging and Parkinson’s Disease" (Drori, Berman, and Mezer, 2021, submitted).
This code automatically computes the main axes of the striatum and other subcortical structures in MRI images, based on the ROI voxel coordinates, using SVD on the individual subject level. It then calculates the spatial functions of MRI intensities (or qMRI values) along these axes on the individual and the group levels. Visualization options are also provided for the ROI axes and the MRI functions.

Use mrGrad_run.m script for an example run.

Requirements:
* MATLAB. Available from: http://www.mathworks.com/products/matlab/
* Vistasoft MATLAB toolbox. Available from: https://github.com/vistalab/vistasoft
Recommended:
* boundedline-pkg MATLAB toolbox. Available from: https://github.com/kakearney/boundedline-pkg