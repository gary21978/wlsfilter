#ifndef _FGS_H_
#define _FGS_H_

int InitFGS(int width, int height);
int DeinitFGS();
int FastGlobalSmoothing(float* image, int width, int height, float sigma, float lambda, 
                        int solver_iteration = 3);
int PrepareBLFKernel(float sigma);

#endif 