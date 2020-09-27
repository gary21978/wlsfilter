#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>

float* buffer1, *buffer2, *buffer3, *N;

#ifdef MATLAB_MEX_FILE  // if compiled as a MEX-file
#include "mex.h"
#include <time.h>
#include "GuidedFilter.h"

// Guided filtering MEX-gateway
// output_image = GuidedFilter(input_image, radius, epsilon);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs < 3) 
    { 
        mexErrMsgTxt("GF must be called with 3 arguments.");
    }
  
    const mxArray *img_input = prhs[0];                          // input images
    mxArray *img_output = plhs[0] = mxDuplicateArray(img_input); // output images
    void *ptr_input, *ptr_output;
    mxClassID category = mxGetClassID(img_input);  
    switch (category)
    {
        case mxUINT8_CLASS: 
            ptr_input = (uint8_t*)mxGetPr(img_input); 
            ptr_output = (uint8_t*)mxGetPr(img_output); break;
        case mxUINT16_CLASS: 
            ptr_input = (uint16_t*)mxGetPr(img_input); 
            ptr_output = (uint16_t*)mxGetPr(img_output); break;    
        case mxSINGLE_CLASS:
            ptr_input = (float*)mxGetPr(img_input);
            ptr_output = (float*)mxGetPr(img_output); break;
        case mxDOUBLE_CLASS: 
            ptr_input = (double*)mxGetPr(img_input);
            ptr_output = (double*)mxGetPr(img_output); break;
        default: mexErrMsgTxt("Expected input image to be one of these types:\n\nuint8, uint16, single, double");
    }

    // image resolution
    int width = mxGetDimensions(img_input)[1];
    int height = mxGetDimensions(img_input)[0];
    int nChannels = (mxGetNumberOfDimensions(img_input) > 2)?(mxGetDimensions(img_input)[2]):1;
    
    // Guided filter parameters
    int radius = mxGetScalar(prhs[1]);
    float epsilon = mxGetScalar(prhs[2]);

    mexPrintf("Image resolution: %d x %d x %d\n", width, height, nChannels);
    mexPrintf("Parameters:\n");
    mexPrintf("    radius = %d\n", radius);
    mexPrintf("    epsilon = %f\n", epsilon);
    
    // Image buffer preperation
    float* image = (float*)malloc(width*height*nChannels*sizeof(float));

    switch (category)
    {
        case mxUINT8_CLASS:
            for (int i = 0; i < height; ++i)
                for (int j = 0; j < width; ++j)
                    for (int c = 0; c < nChannels; ++c)
                        image[c*height*width + i*width + j] = (float)(((uint8_t*)ptr_input)[j*height + i + c*height*width])/255;
            break;
        case mxUINT16_CLASS:
            for (int i = 0; i < height; ++i)
                for (int j = 0; j < width; ++j)
                    for (int c = 0; c < nChannels; ++c)
                        image[c*height*width + i*width + j] = (float)(((uint16_t*)ptr_input)[j*height + i + c*height*width])/65535;
            break;
        case mxSINGLE_CLASS:
            for (int i = 0; i < height; ++i)
                for (int j = 0; j < width; ++j)
                    for (int c = 0; c < nChannels; ++c)
                        image[c*height*width + i*width + j] = (float)(((float*)ptr_input)[j*height + i + c*height*width]);
            break;
        case mxDOUBLE_CLASS:
            for (int i = 0; i < height; ++i)
                for (int j = 0; j < width; ++j)
                    for (int c = 0; c < nChannels; ++c)
                        image[c*height*width + i*width + j] = (float)(((double*)ptr_input)[j*height + i + c*height*width]);
            break;
    }

    for(int i = 0;i < height*width*nChannels; ++i)
    {
        if (image[i]*(image[i] - 1.0) > 0.0)
            mexErrMsgTxt("Intensity value must be between 0 and 1.");
    }
    
    InitGF(width, height, radius);

    clock_t m_begin = clock(); // time measurement;
    for (int c = 0; c < nChannels; ++c)
    {
        GuidedFilter(&image[c*width*height], width, height, radius, epsilon);
    }
    mexPrintf("Elapsed time is %f seconds.\n", double(clock() - m_begin)/CLOCKS_PER_SEC);

    DeinitGF();

    // output conversion
    switch (category)
    {
        case mxUINT8_CLASS:
            for (int i = 0; i < height; ++i)
                for (int j = 0; j < width; ++j)
                    for (int c = 0; c < nChannels; ++c)
                        ((uint8_t*)ptr_output)[j*height + i + c*height*width] = (uint8_t)(image[c*height*width + i*width + j]*255 + 0.5);
            break;
        case mxUINT16_CLASS:
            for (int i = 0; i < height; ++i)
                for (int j = 0; j < width; ++j)
                    for (int c = 0; c < nChannels; ++c)
                        ((uint16_t*)ptr_output)[j*height + i + c*height*width] = (uint16_t)(image[c*height*width + i*width + j]*65535 + 0.5);
            break;
        case mxSINGLE_CLASS:
            for (int i = 0; i < height; ++i)
                for (int j = 0; j < width; ++j)
                    for (int c = 0; c < nChannels; ++c)
                        ((float*)ptr_output)[j*height + i + c*height*width] = (float)(image[c*height*width + i*width + j]);
            break;
        case mxDOUBLE_CLASS:
            for (int i = 0; i < height; ++i)
                for (int j = 0; j < width; ++j)
                    for (int c = 0; c < nChannels; ++c)
                        ((double*)ptr_output)[j*height + i + c*height*width] = (double)(image[c*height*width + i*width + j]);    
            break;    
    }
    free(image);
}
#endif // MATLAB_MEX_FILE

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

/*
 *    Guided filtering Routine
 *    
 *    image                         image buffer of size width*height (float)
 *    width                         image width
 *    height                        image height
 *    radius                        neighborbood radius
 *    epsilon                       variance threshold
*/

int GuidedFilter(float* image, int width, int height, int radius, float epsilon)
{
    int px, sz = width*height;

    if (image == NULL)
    {
        return 1;
    }

    for (px = 0; px < sz; ++px)
    {
        // restrict the pixel value within 0 and 1
        if (image[px]*(image[px] - 1.0) > 0.0)
            return 1;
    }

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

    // clip the value
    for (px = 0; px < sz; ++px)
    {
        if(image[px] < 0.0)
            image[px] = 0.0;
        if(image[px] > 1.0)
            image[px] = 1.0;
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