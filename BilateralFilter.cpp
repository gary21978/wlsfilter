#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>

#ifdef MATLAB_MEX_FILE  // if compiled as a MEX-file
#include "mex.h"
#include <time.h>
#include "BilateralFilter.h"

// Bilateral filtering MEX-gateway
// output_image = BilateralFilter(input_image, sigma_s, sigma_r);

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

    // Bilateral filter parameters
    float sigma_s = mxGetScalar(prhs[1]);
    float sigma_r = mxGetScalar(prhs[2]);

    mexPrintf("Image resolution: %d x %d x %d\n", width, height, nChannels);
    mexPrintf("Parameters:\n");
    mexPrintf("    spatial sigma = %f\n", sigma_s);
    mexPrintf("    range sigma = %f\n", sigma_r);
    
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


    clock_t m_begin = clock(); // time measurement;
    for (int c = 0; c < nChannels; ++c)
    {
        BilateralFilter(&image[c*width*height], width, height, sigma_s, sigma_r);
    }
    mexPrintf("Elapsed time is %f seconds.\n", double(clock() - m_begin)/CLOCKS_PER_SEC);


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

/*
 *    Bilateral filtering
 *    
 *    image                         image buffer of size width*height (float)
 *    width                         image width
 *    height                        image height
 *    sigma_s                       spatial sigma
 *    sigma_r                       range sigma
*/

int BilateralFilter(float* image, int width, int height, float sigma_s, float sigma_r)
{
    int px, sz = width*height;
    int i, j, k, l, imin, imax, jmin, jmax, r;
    float sum1, sum2, bf;

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

    r = ceil(2.0*sigma_s) + 1;
    if (r >= 10)
        r = 10;

    float* spatial_filter = (float*)malloc((2*r + 1)*(2*r + 1)*sizeof(float));
    float* image_bf = (float*)malloc(sz*sizeof(float));

    for (k = -r;k <= r; ++k)
    {
        for (l = -r;l <= r; ++l)
        {
            spatial_filter[(k + r)*(2*r + 1) + l + r] 
            = exp(-(float)(k*k + l*l)/(2.0*sigma_s*sigma_s));
        }
    }

    for (i = 0;i < height; ++i)
    {
        for (j = 0;j < width; ++j)
        {
            imin = (i - r > 0)?(i - r):0;
            imax = (i + r < height)?(i + r):(height - 1);
            jmin = (j - r > 0)?(j - r):0;
            jmax = (j + r < width)?(j + r):(width - 1);

            sum1 = 0.0;
            sum2 = 0.0;
            for (k = imin - i;k <= imax - i; ++k)
            {
                for (l = jmin - j;l <= jmax - j; ++l)
                {
                    bf = image[(i + k)*width + j + l] - image[i*width + j];
                    bf = bf/sigma_r;
                    bf = exp(-0.5*bf*bf);
                    bf *= spatial_filter[(k + r)*(2*r + 1) + l + r];
                    sum1 += image[(i + k)*width + j + l]*bf;
                    sum2 += bf;
                }
            }
            image_bf[i*width + j] = sum1/sum2;
        }
    }
    memcpy(image, image_bf, sz*sizeof(float));

    // clip the value
    for (px = 0; px < sz; ++px)
    {
        if(image[px] < 0.0)
            image[px] = 0.0;
        if(image[px] > 1.0)
            image[px] = 1.0;
    }

    free(spatial_filter);
    free(image_bf);
    return 0;
}