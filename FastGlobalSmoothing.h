void InitFGS(int width, int height);
void DeinitFGS();
void FastGlobalSmoothing(float* image, int width, int height, float sigma, float lambda, 
                         int solver_iteration = 3);
void PrepareBLFKernel(float sigma);
