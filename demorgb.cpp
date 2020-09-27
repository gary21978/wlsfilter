#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image.h"
#include "stb/stb_image_write.h"
#include "FastGlobalSmoothing.h"

int main(int argc, char* argv[])
{
    int width, height, bpp;
    if (argc == 1)
    {
        printf("Usage: %s input [output] [sigma] [lambda]\n",argv[0]);
        return -1;
    }
    uint8_t* image = stbi_load(argv[1], &width, &height, &bpp, 0);
    
    if (image == NULL)
    {
        printf("Cannot open file %s!\n",argv[1]);
        return -1;        
    }
    const int sz = width*height;
    float* channel = (float*)malloc(sizeof(float)*sz*bpp);

    for (int px = 0;px < sz; ++px)
    {
        for (int c = 0; c < bpp; ++c)
        {
            channel[c*sz + px] = (float)image[bpp*px + c];
        }
    }

    const float sigma = (argc > 3)?atof(argv[3]):0.08;
    const float lambda = (argc > 4)?atof(argv[4]):10.0;

    InitFGS(width, height);
    clock_t start = clock();
    for (int c = 0; c < bpp; ++c)
    {
        FastGlobalSmoothing(&channel[c*sz], width, height, sigma, lambda);
    }
    printf("exec time: %.2f ms\n",1000*(clock() - start)/(float)CLOCKS_PER_SEC);
    DeinitFGS();

    for (int px = 0;px < sz; ++px)
    {
        for (int c = 0; c < bpp; ++c)
        {
            image[bpp*px + c] = (u_int8_t)channel[c*sz + px];
        }
    }

    bool correct_format = false;
    if (argc > 2)
    {
        int len = strlen(argv[2]);
        if (len > 4)
        {
            char extension[5];
            strcpy(extension, &argv[2][len - 4]);
            if (strcmp(extension, ".bmp") == 0)
            {
                correct_format = true;
                stbi_write_bmp(argv[2], width, height, bpp, image);
            }
            else 
            {
                if (strcmp(extension, ".png") == 0)
                {
                    correct_format = true;
                    stbi_write_png(argv[2], width, height, bpp, image, width*bpp);
                }
                else
                {
                    if (strcmp(extension, ".jpg") == 0)
                    {   
                        correct_format = true;
                        stbi_write_jpg(argv[2], width, height, bpp, image, 100);
                    }
                }
            }
        }
    }

    if (!correct_format)
    {
        printf("invalid output filename!\n");
        stbi_write_png("output.png", width, height, bpp, image, width*bpp);
    }

    free(channel);
    free(image);
    return 0;
}
