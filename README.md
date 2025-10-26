# **mrGrad: Deep Brain MRI Gradient Analysis Toolbox**

`mrGrad` is a MATLAB-based toolbox for analyzing spatial gradients within deep brain structures using MRI data.
It provides tools to automatically compute principal anatomical axes within subcortical regions and quantify MRI-derived microstructural profiles along those axes.

Originally developed to characterize the **striatal gradient organization** in aging and Parkinson’s disease (*Drori, Berman & Mezer, Science Advances 2022*),
the framework is general and can be applied to **other subcortical structures**.

---

## **Overview**

The toolbox:

* Computes **principal axes** of a region of interest (ROI) using **singular value decomposition (SVD)** on the voxel coordinates of the ROI mask.
* Segments the ROI into **equidistant bins** along the computed axes.
* Calculates **spatial functions (gradients)** of MRI or quantitative MRI (qMRI) values along these axes, both at the individual and group levels.
* Produces **subject-level NIfTI maps** of the axis-based segmentations and **summary outputs** for further statistical analysis.

---

## **Outputs**

`mrGrad` produces:

* ROI **axes definitions** for each subject
* **Segmented NIfTI volumes** aligned to those axes
* **Tables and CSV files** summarizing voxel intensity profiles
* **Group-level summaries** for comparing spatial profiles across subjects or conditions

---

## **Visualization**

Built-in tools allow visualization of:

* ROI axes overlaid on anatomical MRI data
* Spatial profiles (gradients) with optional confidence bounds

*(Plotting with bounded lines can be enhanced using the optional `boundedline` MATLAB package.)*

---

## **Example Usage**

To reproduce the analysis for the striatum, run `mrGrad_run.m`

Additional example scripts demonstrate how the same framework can be applied to other brain regions—
for example, the midbrain, which was developed as an extension of mrGrad in
*Berman, Drori & Mezer, NeuroImage, 2022.*

The midbrain application includes additional code specific to that analysis,
reflecting processing steps specific to this region.

---

## **Requirements**

`mrGrad` is fully self-contained and requires only MATLAB.

**Required:**

* [MATLAB](https://www.mathworks.com/products/matlab/)

**Recommended (for visualization):**

* [`boundedline-pkg`](https://github.com/kakearney/boundedline-pkg)

---

## **Citation**

If you use `mrGrad` in your work, please cite:

**Drori, Berman & Mezer (2022).**
*Mapping Microstructural Gradients of the Human Striatum in Normal Aging and Parkinson’s Disease.*
*Science Advances.*


