Project Title: Hyper-spectral Image Processing

Description:
This project deals with processing hyper-spectral images. Hyper-spectral imaging is a technique where the amount of light reflected from a scene is measured over a wide range of frequencies or wavelengths. In a typical color camera, only three color channels are captured, whereas in hyper-spectral imaging, a large number of color channels are captured. 

When we capture an image, the light that hits various objects in the world enters the camera and is detected by a set of sensors that together represent a pixel. For example, in everyday cameras, light is captured by three sensors, red, green, and blue. Each material reflects a different amount of light at different wavelengths, and when several materials are detected in the same pixel, the total amount of light detected by the sensor is the weighted sum of their contributions. 

The goal of this project is to find the proportion of each material in each pixel of a hyper-spectral image. This can be achieved by solving an optimization problem, where we want to find the material proportions that best match the observed pixel values:

x∗ =argmin ∥y−Hx∥2 x∈Rn
s.t ⃗1Tx=1 x≥0

Here, y is the observed pixel values, H is the matrix that maps material proportions to observed pixel values, and x is the material proportions that we want to find. The optimization problem finds the x that minimizes the difference between the observed pixel values and the pixel values predicted by the material proportions.

Installation:
To use this project, you need Python 3 installed on your computer. You also need the numpy and scipy libraries. These can be installed using pip:

```
pip install numpy
pip install scipy
```

Usage:
The main script for this project is `hyper_spectral_processing.py`. This script takes as input a hyper-spectral image and the corresponding material-to-pixel mapping matrix, and outputs an image where each pixel is color-coded according to the proportions of different materials in that pixel. To run the script, navigate to the directory where the script is located and run:

```
python hyper_spectral_processing.py -i input_image.npy -m material_mapping.npy -o output_image.npy
```

Here, `input_image.npy` is the filename of the input hyper-spectral image, `material_mapping.npy` is the filename of the material-to-pixel mapping matrix, and `output_image.npy` is the filename of the output image.

Contributing:
If you want to contribute to this project, feel free to fork the repository and submit a pull request. We welcome contributions of all kinds, including bug reports, feature requests, and code improvements.

License:
This project is licensed under the MIT License - see the LICENSE file for details.
