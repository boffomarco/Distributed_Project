clear all;
clc;

C = imread("Density100.png");
imshow(C)

BW = rgb2gray(C)
figure,imshow(BW)