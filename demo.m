clear;close all;
img = imread('stripes.png');
filtered = FastGlobalSmoothing(img, 0.03, 900);
subplot(121);imshow(img);title('Original');
subplot(122);imshow(filtered);title('Filtered');