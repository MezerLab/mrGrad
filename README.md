# **mrGrad: Deep Brain MRI Gradient Analysis Toolbox**

`mrGrad` is a MATLAB-based toolbox for analyzing spatial gradients within deep brain structures using MRI data.
It provides tools to automatically compute principal anatomical axes within subcortical regions and quantify MRI-derived microstructural profiles along those axes.

The method was developed for the study:

> **Mapping Microstructural Gradients of the Human Striatum in Normal Aging and Parkinson’s Disease**
> *Drori, Berman, and Mezer — Science Advances, 2022*

A related extension applying the same framework to the **midbrain** was later described in:

> *Berman, Drori, and Mezer — NeuroImage, 2022*
> (implemented in this toolbox as a dedicated example with region-specific scripts)


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


## Midbrain Extension

An extended implementation of the same framework for the **midbrain**
(including specialized preprocessing and example data) is available in:

```
examples/midbrain_extension/
```

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

> Drori, Berman, and Mezer.
> *Mapping Microstructural Gradients of the Human Striatum in Normal Aging and Parkinson’s Disease.*
> *Science Advances, 2022.*

and, if using the midbrain extension:

> Berman, Drori, and Mezer.
> *Spatial profiles provide sensitive MRI measures of the midbrain micro- and macrostructure.*
> *NeuroImage, 2022.*	

