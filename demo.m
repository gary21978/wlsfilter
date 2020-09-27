clear;close all;
img = imread('stripes.png');             % uint8
%img = im2uint16(imread('stripes.png')); % uint16
%img = im2double(imread('stripes.png')); % double
%img = im2single(imread('stripes.png')); % single

filtered = FastGlobalSmoothing(img, 0.03, 900);
%filtered = GuidedFilter(img, 3, 0.25);
%subplot(121);imshow(img);title('Original');
%subplot(122);imshow(filtered);title('Filtered');