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

## Solution approach
The problem at hand is a dual problem, and our chosen approach is to find a solution to the dual problem by minimizing the Lagrangian function. We will be using the Newton method with gradient descent, as the problem is convex, and the constraints are known, resulting in a tight bound on optimization. 

To begin, we developed the Lagrangian function analytically, with the optimal value of the Lagrangian function defined as the function g. Our next step was to determine the λ and μ values that maximize g, while ensuring that λ ≥ 0. We accomplished this using the Newton method with gradient descent, along with the Log Barrier method to handle the constraint λ ≥ 0.

The H^-1 matrix posed a challenge to compute for large dimensions, so we opted to handle the primal problem using the Log Barrier method on the original constraints and gradient descent only. The primal problem involved minimizing the Lagrangian function using the Newton method while handling the constraints with the Log Barrier method.

Overall, this approach allows us to effectively handle the constraints and optimize the problem, resulting in a successful solution to the dual problem.

## Conclusion
This project demonstrates the application of optimization techniques to analyze hyper-spectral images and determine the proportion of materials present in each pixel. The results can have important applications in fields such as remote sensing, environmental monitoring, and medical imaging.
