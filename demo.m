clear;close all;
img = imread('stripes.png');
filtered = FastGlobalSmoothing(double(img), 0.03, 900);
filtered = cast(filtered, class(img));
subplot(121);imshow(img);title('Original');
subplot(122);imshow(filtered);title('Filtered');
