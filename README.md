# Hyper-spectral Camera Pixel Material Analysis

## Overview
This project aims to analyze the materials present in each pixel of a hyper-spectral camera image. In hyper-spectral cameras, each pixel is represented by a vector of intensities captured by multiple sensors corresponding to different materials present in the scene. The goal of this project is to find the percentage of each material present in a given pixel.

## Methodology
Given a pixel vector `y ∈ Rm` and a transition matrix `H ∈ Rm×n` representing the intensities of `n` materials, the objective is to find the proportion of each material in the pixel vector.

The following optimization problem is solved to find the proportion vector `x`:
```
x∗ =argmin ∥y−Hx∥2 x∈Rn
s.t ⃗1Tx=1 x≥0
```
where `x` represents the proportion vector and `⃗1` is a vector of all ones. 

## Code
The code for this project is written in Python. The `main.py` file contains the implementation of the optimization problem. To run the code, first install the necessary packages listed in `requirements.txt`. Then, run `main.py` with the appropriate arguments.

## Conclusion
This project demonstrates the application of optimization techniques to analyze hyper-spectral images and determine the proportion of materials present in each pixel. The results can have important applications in fields such as remote sensing, environmental monitoring, and medical imaging.
