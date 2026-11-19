---
title: "POGO: An R-Based Software Workflow for Automated Leg Stiffness Estimation from C3D Biomechanical Data"

tags:
  - R
  - biomechanics
  - leg stiffness
  - motion capture
  - C3D
  - running
  - countermovement jump

authors:
  - name: Yuming Wang
    affiliation: "1"

  - name: Qichang Mei
    affiliation: "1, 2"
    corresponding: true

affiliations:
  - name: Faculty of Sports Science, Ningbo University, Ningbo, China
    index: 1

  - name: Auckland Bioengineering Institute, The University of Auckland, Auckland, New Zealand
    index: 2

date: 12 November 2026

bibliography: reference.bib
---


# Summary

This work presents POGO, an R-based automated workflow for leg stiffness estimation from synchronized motion-capture and force-platform data stored in C3D files. The workflow reduces the need for manual C3D file processing and biomechanical event identification.

The workflow enables automated C3D data processing, including data extraction, stance-phase detection, biomechanical event identification, force–displacement parameter extraction, leg stiffness calculation, and result export. The software supports batch processing of running and countermovement jump (CMJ) landing biomechanics datasets and is designed to improve the efficiency and reproducibility of leg stiffness analysis.

The workflow was developed for running and CMJ landing studies using three-dimensional motion-capture and force-platform measurements. Demonstrations using running and CMJ landing datasets show how repeated trials can be processed and standardized leg stiffness outputs generated.

# Statement of need

Leg stiffness is an important biomechanical parameter for characterizing lower-limb mechanical behavior during running and CMJ landing. It reflects the ability of the lower-extremity system to resist deformation under external loading and is commonly estimated using spring–mass models [@mcmahon1990].

Traditional leg stiffness estimation relies on synchronized motion-capture and force-platform data. The calculation requires identification of task-specific biomechanical events and extraction of the corresponding force and leg-compression variables during stance [@farley1996]. Different calculation approaches may also influence stiffness estimates depending on the selected kinetic and kinematic variables [@coleman2012].

However, conventional processing workflows often require manual C3D file inspection, event identification, and parameter extraction. These procedures become time-consuming when processing datasets with multiple participants and repeated trials.

POGO integrates C3D data processing, biomechanical event detection, leg stiffness calculation, and standardized result export [@pogo]. It is intended for biomechanics researchers who work with synchronized C3D motion-capture and force-platform recordings and need to process repeated running or CMJ landing trials.

# State of the field

Several software platforms are available for processing biomechanical data from motion-capture and force-platform systems, including commercial and open-source tools for kinematic and kinetic analysis. These platforms provide flexible data-processing capabilities; however, leg stiffness estimation may still require researchers to develop task-specific workflows for event identification and variable extraction.

Existing studies have proposed different approaches to calculate leg stiffness during running and other movement tasks, but these methods are often implemented through customized scripts or manual processing procedures [@coleman2012; @liew2017]. Such approaches may limit reproducibility and increase the workload associated with large biomechanical datasets.

The R package `c3dr` provides functionality for reading and writing C3D motion-capture files [@nolte2026]. POGO builds on C3D file access by integrating force-based contact detection, task-specific event identification, force–displacement parameter extraction, and leg stiffness calculation. It is therefore intended to complement C3D data-reading tools by providing a workflow for subsequent leg stiffness analysis.

# Software design

POGO was developed in R to automate leg stiffness estimation from synchronized three-dimensional motion-capture and force-platform data. The workflow integrates C3D file processing, kinematic and kinetic data extraction, biomechanical event detection, force–displacement parameter extraction, leg stiffness calculation, and standardized result export within a unified framework (see \autoref{fig:workflow}). C3D files are imported and parsed in R, allowing synchronized marker trajectories and ground reaction force (GRF) data to be extracted for subsequent analysis [@nolte2026].

The workflow uses a shared data-import and force-based event-detection stage, followed by task-specific event identification. Force samples are mapped to marker frames using the sampling-rate relationship. Stance or contact is identified using a predefined vertical GRF threshold. Separate event criteria are then applied for running and CMJ landing, allowing the event definition to reflect the movement task rather than imposing one event rule on both tasks. This design reuses the common C3D import and synchronization steps while keeping event rules explicit for each task; extending the workflow to another movement task would require defining and validating an additional event-detection rule.


For running, the lowest vertical position of the right greater trochanter marker is described as the maximum-compression event. The corresponding vertical GRF and marker displacement from initial contact are then used to calculate leg stiffness according to the force–displacement relationship in \autoref{eq:leg-stiffness} [@coleman2012; @farley1996; @mcmahon1990].

$$K_{leg} = \frac{F_{\max}}{\mathrm{\Delta}L}$$


Here, $F_{\max}$ represents the vertical GRF corresponding to the selected event, and $\Delta L$ represents the corresponding marker displacement. The workflow outputs standardized CSV files containing detected event frames, force values, displacement values, stiffness estimates, and processing status, enabling batch analysis of repeated biomechanical trials.


![Workflow of the automated leg stiffness analysis framework. The workflow includes C3D import, kinematic and kinetic data extraction, stance-phase detection, key-event detection, force and displacement extraction, and leg stiffness calculation and result export.\label{fig:workflow}](figures/media/image1.png)

## Using the software

POGO provides a unified workflow for automated leg stiffness analysis across different movement tasks. To demonstrate its application, two previously collected biomechanical datasets comprising running and CMJ landing trials were processed. A common processing framework was applied to both datasets, while task-specific criteria were used to identify the biomechanical events required for the calculation of leg stiffness.

## Illustrative application: Running and CMJ data processing

Two previously collected biomechanical datasets, comprising running and CMJ landing trials, were used to demonstrate the workflow [@wan2026; @wang2026]. C3D files were processed using the common automated framework shown in Code 1. Task-specific biomechanical events were then identified for each movement (see \autoref{fig:events}).

**Code 1. Common force-based event detection**

```r
# POGO: C3D data extraction
d <- pogo_read_c3d(c3d_file)
marker_data <- c3dr::c3d_data(d, format = "wide")
fp <- d$forceplatform[[forceplate_index]]
force_matrix <- as.matrix(fp$forces)
force_raw <- as.numeric(force_matrix[, 3])
force_abs <- abs(force_raw)

# POGO: Common force-based event detection
peak_force_index <- which.max(
  ifelse(is.finite(force_abs), force_abs, -Inf)
)
force_present <- is.finite(force_abs) &
  force_abs > force_threshold
```

### Running

For running, the stance phase is identified from the vertical GRF signal. The lowest vertical position of the right greater trochanter marker is defined as the maximum-compression event. The corresponding force and marker displacement are extracted for leg stiffness calculation (Code 2).

**Code 2. Running: Maximum-compression event detection**

```r
# POGO: Running maximum compression event detection
onset_force_index <- peak_force_index
while (onset_force_index > 1L &&
       isTRUE(force_present[onset_force_index - 1L])) {
  onset_force_index <- onset_force_index - 1L
}

offset_force_index <- peak_force_index
while (offset_force_index < length(force_present) &&
       isTRUE(force_present[offset_force_index + 1L])) {
  offset_force_index <- offset_force_index + 1L
}

onset_row <- .pogo_force_to_marker_row(
  onset_force_index, analog_per_frame, n_marker_frames
)
offset_row <- .pogo_force_to_marker_row(
  offset_force_index, analog_per_frame, n_marker_frames
)

z_col <- marker_cols[[3]]
z_onset <- as.numeric(marker_data[[z_col]][onset_row])
contact_rows <- seq.int(onset_row, offset_row)
z_contact <- as.numeric(marker_data[[z_col]][contact_rows])

displacement <- z_contact - z_onset
abs_displacement <- abs(displacement)
abs_displacement[!is.finite(abs_displacement)] <- -Inf
local_index <- which.max(abs_displacement)
compression_row <- contact_rows[[local_index]]

z_compression <- as.numeric(marker_data[[z_col]][compression_row])
z_displacement <- z_compression - z_onset
```

### CMJ landing

For CMJ landing, the contact phase is identified from the force signal. Peak vertical GRF is used as the reference event, and the corresponding right greater trochanter displacement is extracted from the synchronized kinematic data (Code 3).

**Code 3. CMJ: Peak-force-associated displacement detection**

```r
# POGO: CMJ peak-force-associated displacement detection
onset_force_index <- peak_force_index
while (onset_force_index > 1L &&
       isTRUE(force_present[onset_force_index - 1L])) {
  onset_force_index <- onset_force_index - 1L
}

onset_row <- .pogo_force_to_marker_row(
  onset_force_index, analog_per_frame, n_marker_frames
)
peak_row <- .pogo_force_to_marker_row(
  peak_force_index, analog_per_frame, n_marker_frames
)

z_col <- marker_cols[[3]]
z_onset <- as.numeric(marker_data[[z_col]][onset_row])
compression_row <- peak_row
z_compression <- as.numeric(marker_data[[z_col]][compression_row])
z_displacement <- z_compression - z_onset
```

The workflow generates standardized trial-level outputs for both tasks, including event frames, force values, displacement values, and leg stiffness estimates. These outputs are available for subsequent statistical analysis.

![Representative biomechanical data acquisition and task-specific event identification for running and CMJ landing.\label{fig:events}](figures/media/image2.png)

# Research impact statement

The workflow provides a practical and reproducible approach for quantifying leg stiffness from synchronized motion-capture and force-platform data. By automating biomechanical event detection and the extraction of force and displacement variables, the software reduces the need for time-consuming manual processing and facilitates the analysis of biomechanical data.

The workflow incorporates task-specific event-detection procedures for running and CMJ landing, allowing leg stiffness to be calculated using event definitions appropriate for each movement task. Its batch-processing capability supports the analysis of multiple trials across these tasks, providing a consistent approach to leg stiffness assessment in biomechanical research.


# AI usage disclosure

No AI tools were used during the development of this software.


# Acknowledgements

This work was supported by the National Natural Science Foundation of China (12202216 and 12672425) and the “Mechanics+” Interdisciplinary Top Innovative Youth Fund Project of Ningbo University (GC2024006).

# References


