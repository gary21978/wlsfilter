# Fast Global Smoothing Based on Weighted Least-Squares

Reimplementation of the paper "D. Min, S. Choi, J. Lu, B. Ham, K. Sohn, and M. N. Do, Fast Global Image Smoothing Based on Weighted Least Squares, *IEEE Trans. on Image Processing*, 23(12), 5638-5653, 2014"



![filtering result](./result.png)



## Usage

Two demo codes in MATLAB and C are provided.

### MATLAB API


output_image = FastGlobalSmoothing(input_image, sigma, lambda)

- The input image can be one of these types: **uint8**, **uint16**, **single** or **double**.
- The output image has the same size and data type as the input image.
- If sigma value is negative or zero, then an *adaptive strategy* based on local noise variance estimation will be adopted.
- The binary MEX files for the operation system Linux 64-bit and Windows 64-bit are provided, with extensions mexa64 and mexw64 respectively.

### C API

int FastGlobalSmoothing(float* image, int width, int height, float sigma, float lambda, int solver_iteration = 3)

- The elements of **single channel** data *image* are arranged in **row-major** order.

- The output image buffer **overwrites** the input image buffer.
- The value of input image is assumed to be in the range [0, 1].
  

## License

Copyright (c) 2020, Li Chen All rights reserved.

For research and education purpose only.
