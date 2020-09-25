#include <stdlib.h>
#include <string.h>
#define NORMALIZE_CONSTANT 0.20888568955258335346769626994501

int BoxFilter(float* dst, const float* src, int width, int height, int r)
{
    int sz = width*height;
    float* cumsum = (float*)malloc(sz*sizeof(float));
    int i, j;
    memcpy(cumsum, src, sz*sizeof(float));
    for(i = 1; i < height; ++i)
    {
        for(j = 0; j < width; ++j)
        {
            cumsum[i*width+j] += cumsum[(i - 1)*width + j];
        }
    }

    for(j = 0; j < width; ++j)
    {
        for (i = 0; i <= r; ++i)
        {
            dst[i*width + j] = cumsum[(i + r)*width + j];
        }
        for (i = r + 1; i <= height - r - 1; ++i)
        {
            dst[i*width + j] = cumsum[(i + r)*width + j] - cumsum[(i - r - 1)*width + j];
        }
        for (i = height - r; i <= height - 1; ++i)
        {
            dst[i*width + j] = cumsum[(height - 1)*width + j] - cumsum[(i - r - 1)*width + j];
        }
    }

    memcpy(cumsum, dst, sz*sizeof(float));
    for(i = 0; i < height; ++i)
    {
        for(j = 1; j < width; ++j)
        {
            cumsum[i*width+j] += cumsum[i*width + j - 1];
        }
    }

    for(i = 0; i < height; ++i)
    {
        for (j = 0; j <= r; ++j)
        {
            dst[i*width + j] = cumsum[i*width + j + r];
        }
        for (j = r + 1; j <= width - r - 1; ++j)
        {
            dst[i*width + j] = cumsum[i*width + j + r] - cumsum[i*width + j - r - 1];
        }
        for (j = width - r; j <= width - 1; ++j)
        {
            dst[i*width + j] = cumsum[i*width + width - 1] - cumsum[i*width + j - r - 1];
        }
    }
    free(cumsum);
    return 0;
}

float NoiseVarianceEstimation(float* variance, const float* img, int width, int height, int r)
{
    float * convolution = (float*)malloc(width*height*sizeof(float));
    float value;
    
    memset(convolution, 0, width*height*sizeof(float));

    //                      [ 1,-2, 1]
    // convolve with kernel [-2, 4,-2]
    //                      [ 1,-2, 1]

    for (int i = 1; i < height - 1; ++i)
    {
        for (int j = 1; j < width - 1; ++j)
        {
            value = 4*img[i*width + j];
            value += img[(i - 1)*width + j - 1];
            value += img[(i - 1)*width + j + 1];
            value += img[(i + 1)*width + j - 1];
            value += img[(i + 1)*width + j + 1];
            value -= 2*img[i*width + j - 1];
            value -= 2*img[i*width + j + 1];
            value -= 2*img[(i - 1)*width + j];
            value -= 2*img[(i + 1)*width + j];
            if (value < 0)
                value *= -1;
            convolution[i*width + j] = value;
        }
    }

    BoxFilter(variance, convolution, width, height, r);

    int w2, h2;
    float total_variance = 0.0;
    for (int i = 0; i < height; ++i)
    {
        h2 = (i < r)?(r - 1 + i):((i >= height - r)?(height + r - i - 2):(2*r - 1));
        for (int j = 0; j < width; ++j)
        {
            w2 = (j < r)?(r - 1 + j):((j >= width - r)?(width + r - i - 2):(2*r - 1));
            variance[i*width + j] *= NORMALIZE_CONSTANT/(double)w2/(double)h2;
            variance[i*width + j] *= variance[i*width + j];
            total_variance += convolution[i*width + j];
        }
    }
    total_variance *= NORMALIZE_CONSTANT/(double)(width - 2)/(double)(height - 2);
    total_variance *= total_variance;

    free(convolution);
    return total_variance;
}