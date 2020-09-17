#include <math.h>
#include <string.h>

#define MAX_VAL 255
float BLFKernel[MAX_VAL + 1];
float *a1_vec, *b1_vec, *c1_vec, *a2_vec, *b2_vec, *c2_vec;

#ifdef MATLAB_MEX_FILE  // if compiled as a MEX-file
#include "mex.h"
#include <time.h>
#include "FastGlobalSmoothing.h"

// mex call
// output_image = FastGlobalSmoothing(input_image, sigma, lambda, solver_iteration);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs < 3) 
    { 
        mexErrMsgTxt("FGS must be called with 3 arguments.");
    }
  
    const mxArray *img = prhs[0]; // input images

    // image resolution
    int width = mxGetDimensions(img)[1];
    int height = mxGetDimensions(img)[0];
    int nChannels = (mxGetNumberOfDimensions(img) > 2)?(mxGetDimensions(img)[2]):1;
    
    // FGS parameters
    float sigma = mxGetScalar(prhs[1]);
    float lambda = mxGetScalar(prhs[2]);
    int solver_iteration = (nrhs > 3)?(int)mxGetScalar(prhs[3]):3;

    mexPrintf("Image resolution: %d x %d x %d\n", width, height, nChannels);
    mexPrintf("Parameters:\n");
    mexPrintf("    Sigma = %f\n", sigma);
    mexPrintf("    Lambda = %f\n", lambda);
    mexPrintf("    Iteration = %d\n", solver_iteration);
    
    // Image buffer preperation
    float* image = (float*)malloc(width*height*nChannels*sizeof(float));
    double* ptr_input = (double*)mxGetPr(img);
    for (int i = 0; i < height; ++i)
        for (int j = 0; j < width; ++j)
            for (int c = 0; c < nChannels; ++c)
                image[c*height*width + i*width + j] = (float)(ptr_input[j*height + i + c*height*width]);
    
    InitFGS(width, height);

    clock_t m_begin = clock(); // time measurement;
    for (int c = 0; c < nChannels; ++c)
    {
        FastGlobalSmoothing(&image[c*width*height], width, height, sigma, lambda, solver_iteration);
    }
    mexPrintf("Elapsed time is %f seconds.\n", double(clock() - m_begin)/CLOCKS_PER_SEC);

    DeinitFGS();

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
 *   Prepare lookup table of bilateral filter kernel
 */

int PrepareBLFKernel(float sigma)
{
    for (int i = 0;i <= MAX_VAL;i++)
    {
        float r = i/(float)MAX_VAL/sigma;
        BLFKernel[i] = exp(-0.5*r*r);     
    }
    return 0;
}

/*
 *    Fast Global Smoothing Routine
 *    
 *    image                         image buffer of size width*height (float)
 *    width                         image width
 *    height                        image height
 *    sigma                         range standard deviation
 *    lambda                        regularization parameter
 *    solver_iteration(optional)    number of iterations (default: 3)
*/

int FastGlobalSmoothing(float* image, int width, int height, 
                        float sigma, float lambda, int solver_iteration)
{
    float *pc, *po, *pa_vec, *pb_vec, *pc_vec;
    float tcr, ncr, tpr, weight, m, lambda_in;
    int i, j, iter, range;

    switch (solver_iteration)
    {
        case 1:  lambda_in = 0.5000*lambda;break;
        case 2:  lambda_in = 0.4000*lambda;break;
        case 3:  lambda_in = 0.3810*lambda;break;
        case 4:  lambda_in = 0.3765*lambda;break;
        default: lambda_in = 0.3750*lambda;break;
        //lambda_in = 1.5*lambda*pow(4.0, solver_iteration - 1)/(pow(4.0, solver_iteration) - 1.0);break;
    }

    // prepare LUT
    PrepareBLFKernel(sigma); 

    for(iter = 0;iter < solver_iteration;iter++) 
    {
        // apply 1D solver along rows of image
        for(i = 0;i < height;i++) 
        {
            memset(a1_vec, 0, sizeof(float)*width);
            memset(b1_vec, 0, sizeof(float)*width);
            memset(c1_vec, 0, sizeof(float)*width);

            pc = &image[i*width];
            tpr = *pc++;      
            pa_vec = &a1_vec[1];
            pc_vec = &c1_vec[0];

            // prepare tridiagonal system
            for(j = 1;j < width;j++) 
            { 
                tcr = *pc++;
                range = (tpr > tcr)?(tpr - tcr):(tcr - tpr);
                weight = BLFKernel[range];
                                
                *pa_vec = -lambda_in*weight;
                *pc_vec = *pa_vec; 
                pa_vec++; 
                pc_vec++;
                tpr = tcr;
            }

            pa_vec = a1_vec;
            pc_vec = c1_vec;
            pb_vec = b1_vec;
            for(j = 0;j < width;j++) 
            {
                *pb_vec++ = 1.0 - *pa_vec++ - *pc_vec++;
            }
            
            // tridiagonal solver
            c1_vec[0] /= b1_vec[0];
            image[i*width] /= b1_vec[0];
            
            pb_vec = &b1_vec[1];
            pa_vec = &a1_vec[1];
            pc_vec = &c1_vec[0];

            pc = &image[i*width];
            tpr = *pc++;
            po = &image[i*width + 1];
            for (j = 1; j < width; j++) 
            {
                m = 1.0/(*pb_vec - *pa_vec*(*pc_vec++));
                *pc_vec *= m; 
                tcr = *pc++;
                ncr = (tcr - *pa_vec*tpr)*m;
                pb_vec++; 
                pa_vec++;
                tpr = ncr;
                *po++ = ncr;
            }

            pc_vec = &c1_vec[width - 2];
            tpr = *--pc;
            --po;
            
            for (j = width - 1; j-- > 0; ) 
            {
                tcr = *--pc;
                ncr = tcr - *pc_vec*tpr;
                pc_vec--;
                tpr = ncr;
                *--po = ncr;
            }
        }

        // apply 1D solver along columns
        for(j = 0;j < width;j++) 
        {
            memset(a2_vec, 0, sizeof(float)*height);
            memset(b2_vec, 0, sizeof(float)*height);
            memset(c2_vec, 0, sizeof(float)*height);

            pc = &image[j];
            tpr = *pc++;
            pa_vec = &a2_vec[1];
            pc_vec = &c2_vec[0];

            // prepare tridiagonal system
            for(i = 1;i < height;i++) 
            {
                pc += (width - 1);
                tcr = *pc++;
                range = (tpr > tcr)?(tpr - tcr):(tcr - tpr);
                weight = BLFKernel[range];

                *pa_vec = -lambda_in*weight;                
                *pc_vec = *pa_vec; 
                pa_vec++; 
                pc_vec++;
                tpr = tcr;
            }

            pa_vec = a2_vec;
            pc_vec = c2_vec;
            pb_vec = b2_vec;
            for(i = 0;i < height;i++) 
            {
                *pb_vec++ = 1.0 - *pa_vec++ - *pc_vec++;
            }

            // tridiagonal solver
            c2_vec[0] /= b2_vec[0];
            image[j] /= b2_vec[0];

            pb_vec = &b2_vec[1];
            pa_vec = &a2_vec[1];
            pc_vec = &c2_vec[0];

            pc = &image[j];
            tpr = *pc++;
            po = &image[j + width];
            for (i = 1; i < height; i++)
            {
                m = 1.0/(*pb_vec - *pa_vec*(*pc_vec++));
                *pc_vec *= m;
                pc += (width - 1);
                tcr = *pc++;
                ncr = (tcr - *pa_vec*tpr )*m;
                pb_vec++; 
                pa_vec++;
                tpr = ncr;
                *po = ncr;
                po += width;
            }

            pc_vec = &c2_vec[height - 2];
            tpr = *--pc;            
            po -= (2*width - 1);
                        
            for (i = height - 1; i-- > 0; ) 
            {
                pc -= (width - 1); 
                tcr = *--pc;
                ncr = tcr - *pc_vec*tpr;
                pc_vec--;
                tpr = ncr; 
                *--po = ncr;
                po -= (width - 1);
            }
        }
        lambda_in *= 0.25;
    }
    return 0;
}

int InitFGS(int width, int height)
{
    a1_vec = (float *)malloc(sizeof(float)*width);
    b1_vec = (float *)malloc(sizeof(float)*width);
    c1_vec = (float *)malloc(sizeof(float)*width);  
    a2_vec = (float *)malloc(sizeof(float)*height);
    b2_vec = (float *)malloc(sizeof(float)*height);
    c2_vec = (float *)malloc(sizeof(float)*height);
    if (a1_vec == NULL || b1_vec == NULL || c1_vec == NULL 
     || a2_vec == NULL || b2_vec == NULL || c2_vec == NULL)
    {
        return 1;
    }
    return 0;
}

int DeinitFGS()
{
    free(a1_vec); free(b1_vec); free(c1_vec);
    free(a2_vec); free(b2_vec); free(c2_vec);
    return 0;
}
