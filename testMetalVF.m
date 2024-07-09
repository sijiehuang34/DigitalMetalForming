%% testMetalVF.m
clc; clear; close all;

% Load the image matrix
load("Acquisition500.mat");
inputImgMat = img;

% Turn into grayscale for pixel histogram analysis
gray_img = im2gray(inputImgMat);

% Binarize (retains top side of the punch only)
bw = imbinarize(gray_img); 

% Denoise (removes small objects)
minSize = 175; % pxl threshold through trial & error
bw = bwareaopen(bw,minSize);

% Make known, fixed-position, irrelevant regions = 0 pxl
bw(1:60, :) = 0;
bw(1:200, 500:end) = 0;

% Now denoise again (removes any remaining small objects)
minSize = 175; % pxl threshold through trial & error
bw = bwareaopen(bw,minSize);

% Fill holes (regionprops is used to estimate the area enclosed later)
% DEBUG NOTES: imfill currently doesn't support string arguments
% Must move this out of the entry-point function
bw = imfill(bw,"holes");

% Look for exterior boundaries
    % "noholes" = don't search for inner contours
    % B = a cell storing the boundary coordinates for every bounded shape
    % L = a label matrix of each shape
[B,L] = bwboundaries(bw,"noholes");

% Image Processing Function
metalVisionFeedback(bw, B, L);