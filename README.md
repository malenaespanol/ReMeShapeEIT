# ReMeShapeEIT

This repository contains a MATLAB implementation of the methods developed in the paper

> **Qualitative reconstruction methods for imaging interior Robin interfaces in EIT from Robin-to-Dirichlet data**  
> **Rafael Ceja Ayala**, **Malena I. Español**, **Govanni Granados**

---

## Overview

This repository contains two **MATLAB files**, `ReMEShapeEIT_CircularInterface.m and ReMEShapeEIT_StarInterface.m  which reproduces the numerical experiments from the paper in Sections 4.1 and 4.2 (in the last (published) version), respectively.

The code implements two qualitative (non-iterative) reconstruction methods for Electrical Impedance Tomography (EIT):

- **Linear Sampling Method (LSM)**
- **Regularized Factorization Method (RFM)**

using three different regularization techniques

- **Tikhonov Regularization**
- **Spectral Cut off**
- **Truncated Total Least Squars (TTLS)**

The setting involves Robin boundary conditions on both the exterior boundary and an interior interface, and reconstruction is performed using the Robin-to-Dirichlet (RtD) map. The code generates indicator functions for the interior region and reproduces the figures in the paper, including noise-free and noisy data, different interior region sizes, and various noise levels and boundary sampling resolutions.

---

## Citation

If you use this code, please cite the associated paper:

@article{ceja2026reconstruction,
  title={Qualitative reconstruction methods for imaging interior Robin interfaces in EIT from Robin-to-Dirichlet data},
  author={Ceja Ayala, Rafael and Español, Malena I. and Granados, Govanni},
  journal={Inverse Problems},
  year={2026}
}
