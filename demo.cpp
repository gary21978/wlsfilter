#include <stdlib.h>
#include <stdio.h>
#include "FastGlobalSmoothing.h"

int main(int argc, char* argv[])
{
    const int width = 1920, height = 1080;
    const int sz = width*height;
    u_int8_t* buffer  = (u_int8_t*)malloc(sz*3/2);
    float* image = (float*)malloc(sizeof(float)*sz);
    // Load YUV data from file

    FILE* fin = fopen("input.yuv","rb");
    if (fin == NULL)
    {
        abort();
    }
    int ret = fread(buffer, sz*3/2, 1, fin);
    fclose(fin);

    // convert luma component
    for(int px = 0;px < sz; ++px)
    {
        image[px] = (float)buffer[px];
    }

    const float sigma = 0.05;
    const float lambda = 10.0;

	// fast global smoothing
    InitFGS(width, height);
    FastGlobalSmoothing(image, width, height, sigma, lambda);
    DeinitFGS();
    
    // overwrite luma component
    for(int px = 0;px < sz; ++px)
    {
        buffer[px] = (u_int8_t)(image[px]);
    }

    // write filtered to file
    FILE* fout = fopen("output.yuv","wb");
    if (fout == NULL)
    {
        abort();
    }
    fwrite(buffer, sz*3/2, 1, fout);
    fclose(fout);

    free(buffer);
    free(image);
    return 0;
}
