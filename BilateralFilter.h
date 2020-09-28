#ifndef _BF_H_
#define _BF_H_

// naive bilateral filter

int BilateralFilter(float* image, int width, int height, float sigma_s, float sigma_r);

#endif