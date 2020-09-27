#ifndef _GF_H_
#define _GF_H_

int GuidedFilter(float* image, int width, int height, int radius, float epsilon);
int InitGF(int width, int height, int radius);
int DeinitGF();

#endif 