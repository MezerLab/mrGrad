# mrGrad: Deep Brain MRI Gradient Analysis Toolbox (v2.0)

`mrGrad` is a MATLAB-based toolbox for analyzing spatial gradients within deep brain structures using MRI data.
It automatically computes principal anatomical axes within subcortical regions and quantifies MRI-derived microstructural profiles along those axes.

Developed for:

**Mapping Microstructural Gradients of the Human Striatum in Normal Aging and Parkinson’s Disease**
*Drori, Berman, and Mezer. Science Advances, 2022*

An extended implementation for the midbrain (including additional region-specific code and example data) was later described in:
**Spatial profiles provide sensitive MRI measures of the midbrain micro- and macrostructure**
*Berman, Drori, and Mezer. NeuroImage, 2022*

---

## Overview

The toolbox performs a full pipeline for regional gradient analysis in MRI data:

1. **ROI Axis Extraction**
   Derives principal axes of subcortical ROIs (e.g., caudate, putamen) using singular value decomposition (SVD) on voxel coordinates.
   Axes represent dominant anatomical directions (e.g., anterior–posterior, dorsal–ventral, medial–lateral).

2. **ROI Segmentation**
   Divides each ROI into bins along each axis.
   Supports segmentation methods:

   * `equidistance` — equal spatial spacing (default)
   * `equivolume` — equal voxel volume per bin

3. **Gradient Computation**
   Samples voxel values from one or more MRI parameter maps (e.g., R1, MTsat, Water Fraction) along each axis.
   Computes per-segment statistics (mean or median).
   Supports multiple modalities per subject (v2.0 feature).

4. **Group-Level Summary**
   Aggregates across subjects, producing mean, standard deviation, and standard error of the mean (SEM) profiles.
   Outputs both MATLAB structs and CSV tables for flexible downstream analysis.

---

## Inputs

`mrGrad` expects a MATLAB struct `Data` with subject-level file paths and metadata:

| Field                          | Description                                                                    |
| ------------------------------ | ------------------------------------------------------------------------------ |
| `seg_list`                     | {N×1} cell array of segmentation or mask NIfTI files (one per subject)         |
| `map_list`                     | {N×M} cell array of coregistered MRI parameter maps (M parameters per subject) |
| `subject_ids` *(optional)*     | {N×1} subject identifiers (e.g., `'sub-001'`, `'sub-002'`)                     |
| Additional fields *(optional)* | e.g., `Age`, `Group`, `Sex` — retained in results                              |

Example:

```matlab
Data.seg_list = {'sub1_seg.nii.gz'; 'sub2_seg.nii.gz'; 'sub3_seg.nii.gz'};
Data.map_list = {
  'sub1_R1.nii.gz', 'sub1_MTsat.nii.gz';
  'sub2_R1.nii.gz', 'sub2_MTsat.nii.gz';
  'sub3_R1.nii.gz', 'sub3_MTsat.nii.gz'};
Data.subject_ids = {'sub-1'; 'sub-2'; 'sub-3'};
```

---

## Key Options

| Parameter             | Description                                                  |
| --------------------- | ------------------------------------------------------------ |
| `'ROI'`               | Vector of label indices (e.g., `[11, 12, 50, 51]`)           |
| `'roi_names'`         | ROI names corresponding to labels                            |
| `'Axes'` or `'PC'`    | Anatomical ROI axes to analyze (default: `[1 2 3]`)          |
| `'n_segments'`        | Number of bins per axis (default: 7)                         |
| `'segmenting_method'` | `'equidistance'` (default) or `'equivolume'`                 |
| `'stat'`              | Statistic used per segment: `'median'` (default) or `'mean'` |
| `'parameter_names'`   | Names of MRI parameters (e.g., `{'R1','MTsat'}`)             |
| `'units'`             | Units of each parameter (e.g., `{'1/s','p.u.'}`)             |
| `'output_mode'`       | `'minimal'`, `'default'`, or `'extended'` (see below)        |
| `'Parallel'`          | Use MATLAB Parallel Toolbox (optional)                       |

---

## Outputs

`mrGrad` produces both MATLAB and file-based outputs, depending on the chosen `output_mode`:

### 1. `minimal`

* Saves only summary statistics:

  * `mrGrad_out.mat` — compact summary struct (subject- and group-level gradient results)
  * `mrGrad_out.csv` — combined summary table across all ROIs and subjects
* No subject-level axes data or NIfTI outputs.

### 2. `default`

* Saves:

  * `mrGrad_out.mat` — full results (including subject-level axes data)
  * `mrGrad_out.csv` — summary table suitable for statistical analysis
* Includes subject-level axes data in the MATLAB struct but does not save NIfTI masks.

### 3. `extended`

* Includes everything from *default* mode, plus:

  * 3D NIfTI segmentation files per subject and ROI segment.
    Each corresponds to one axis (e.g., `sub-1/mrGradSeg_left-putamen_axis1.nii.gz`).
* Recommended for visual inspection and post-hoc analyses.

---

### Example Output Structure

```
ExampleResults/
│
├── mrGrad_out.mat
├── mrGrad_out.csv
└── mrGradSeg/         (only in 'extended' mode)
    ├── sub-1/mrGradSeg_left-putamen_axis1.nii.gz
    ├── sub-1/mrGradSeg_left-putamen_axis2.nii.gz
    └── ...
```

Each ROI result (`RG.Left_Putamen`, `RG.Right_Caudate`, etc.) contains:

* `Results`: per-axis, per-parameter gradient profiles (subject- and group-levels)
* `individual_data`: per-subject axes data (if not minimal)
* `user_input_fields`: additional input metadata (e.g., age, group)

---

## Example: Striatum Analysis

Run the included example:

```matlab
mrGrad_run
```

This executes a full analysis for the bilateral caudate and putamen using multiple MRI parameters (`R1`, `WaterFraction`, `MTsat`).

Outputs will be saved to:

```
example_data/ExampleResults/
```

---

## Visualization and Post-Processing

| Function                | Description                                       |
| ----------------------- | ------------------------------------------------- |
| `mrgrad_show_gradients` | Plot individual and group-level gradient curves   |
| `mrgrad_axis_visualize` | Visualize ROI segmentations and axes in MRI space |
| `mrGrad_average_LR`     | Compute mean left–right ROI gradients             |
| `mrGrad_subset`         | Select subsets of subjects                        |
| `mrGrad_asymmetry`      | Compute inter-hemispheric asymmetry indices       |

Example:

```matlab
mrgrad_show_gradients(RG,'error_name','SEM');
RG_avg = mrGrad_average_LR(RG);
mrgrad_show_gradients(RG_avg,'error_name','SEM');
```

---

## Midbrain Extension

An extension code adapted for the midbrain, is available in:

```
workflows/MidBrainProfiles/
```
Run:

```matlab
Calc_MB_profile_example_data.m
```
---

## Requirements

**Required:**

* [MATLAB](https://www.mathworks.com/products/matlab.html)

**Recommended:**

* [`boundedline-pkg`](https://github.com/kakearney/boundedline-pkg) for plotting bounded mean±error curves

---

## Version 2.0 Highlights

* Support for multiple MRI parameters per subject
* Flexible output modes (`minimal`, `default`, `extended`)
* Enhanced visualization utilities and summary exports

---

## Citation

If you use `mrGrad` toolbox or its underlying conceptual framework in your work, please cite:

Drori, Berman, and Mezer.
*Mapping Microstructural Gradients of the Human Striatum in Normal Aging and Parkinson’s Disease.*
*Science Advances, 2022.*

and, if using the midbrain extension:

Berman, Drori, and Mezer.
*Spatial profiles provide sensitive MRI measures of the midbrain micro- and macrostructure.*
*NeuroImage, 2022.*


