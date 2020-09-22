#include <math.h>
#include <string.h>
#include <stdio.h>

float* buffer1, *buffer2, *buffer3, *N;

#ifdef MATLAB_MEX_FILE  // if compiled as a MEX-file
#include "mex.h"
#include <time.h>

// mex call
// output_image = FastGlobalSmoothing(input_image, sigma, lambda, solver_iteration);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs < 3) 
    { 
        mexErrMsgTxt("GF must be called with 3 arguments.");
    }
  
    const mxArray *img = prhs[0]; // input images

    // image resolution
    int width = mxGetDimensions(img)[1];
    int height = mxGetDimensions(img)[0];
    int nChannels = (mxGetNumberOfDimensions(img) > 2)?(mxGetDimensions(img)[2]):1;
    
    // Guided filter parameters
    int radius = mxGetScalar(prhs[1]);
    float epsilon = mxGetScalar(prhs[2]);

    mexPrintf("Image resolution: %d x %d x %d\n", width, height, nChannels);
    mexPrintf("Parameters:\n");
    mexPrintf("    radius = %f\n", radius);
    mexPrintf("    epsilon = %f\n", epsilon);
    
    // Image buffer preperation
    float* image = (float*)malloc(width*height*nChannels*sizeof(float));
    double* ptr_input = (double*)mxGetPr(img);
    for (int i = 0; i < height; ++i)
        for (int j = 0; j < width; ++j)
            for (int c = 0; c < nChannels; ++c)
                image[c*height*width + i*width + j] = (float)(ptr_input[j*height + i + c*height*width]);
    
    InitGF(width, height, radius);

    clock_t m_begin = clock(); // time measurement;
    for (int c = 0; c < nChannels; ++c)
    {
        GuidedFilter(&image[c*width*height], width, height, radius, epsilon);
    }
    mexPrintf("Elapsed time is %f seconds.\n", double(clock() - m_begin)/CLOCKS_PER_SEC);

    DeinitGF();

    // output
    mxArray *image_result = plhs[0] = mxDuplicateArray(img);
    double* ptr_output = (double*)mxGetPr(image_result);

    for (int i = 0; i < height; ++i)
        for (int j = 0; j < width; ++j)
            for (int c = 0; c < nChannels; ++c)
                ptr_output[j*height + i + c*height*width] = (double)(image[c*height*width + i*width + j]);

    free(image);
}
#endif // MATLAB_MEX_FILE

/*
 *    Guided Filter Routine
 *    
 *    image                         image buffer of size width*height (float)
 *    width                         image width
 *    height                        image height
 *    radius                        local window radius
 *    epsilon                       regularization parameter
*/

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

int GuidedFilter(float* image, int width, int height, int radius, float epsilon)
{
    int sz = width*height;
    int px;
    for (px = 0; px < sz; ++px)
    {
        buffer1[px] = image[px]*image[px];
    }
    BoxFilter(buffer2, buffer1, width, height, radius);
    for (px = 0; px < sz; ++px)
    {
        buffer2[px] /= N[px];  // buffer2:mean_II
    }
    BoxFilter(buffer1, image, width, height, radius);
    for (px = 0; px < sz; ++px)
    {
        buffer1[px] /= N[px];  // buffer1:meanI
    }
    for (px = 0; px < sz; ++px)
    {
        buffer2[px] -= buffer1[px]*buffer1[px];  // buffer2: varI
        buffer2[px] /= (buffer2[px] + epsilon);  // buffer2: a
        buffer1[px] *= (1.0 - buffer2[px]);      // buffer1: b
    }

    BoxFilter(buffer3, buffer2, width, height, radius);
    for (px = 0; px < sz; ++px)
    {
        buffer3[px] /= N[px];  // buffer3:meana
    }

    BoxFilter(buffer2, buffer1, width, height, radius);
    for (px = 0; px < sz; ++px)
    {
        buffer2[px] /= N[px];  // buffer2:meanb
    }
    
    for (px = 0; px < sz; ++px)
    {
        image[px] = buffer3[px]*image[px] + buffer2[px];
    }
    return 0;
}

int InitGF(int width, int height, int radius)
{
    int sz = width*height;
    N = (float*)malloc(sz*sizeof(float));
    buffer1 = (float*)malloc(sz*sizeof(float));
    buffer2 = (float*)malloc(sz*sizeof(float));
    buffer3 = (float*)malloc(sz*sizeof(float));
    // size of each local patch
    float* ones = (float*)malloc(sz*sizeof(float));
    for (int px = 0;px < sz;++px)
    {
        ones[px] = 1.0;
    }
    BoxFilter(N, ones, width, height, radius);
    free(ones);
    return 0;
}

int DeinitGF()
{
    free(buffer1);
    free(buffer2);
    free(buffer3);
    free(N);
}