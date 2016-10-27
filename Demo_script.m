% Demo_script.

%% Clear everyting existing.

clear ;
close all;
clc;

%% Load Data.

load('GradientData.mat','x','y','z','tsx','tsy','nsx','nsy');

% Choose the slope to integrate.
sx = tsx;
sy = tsy;

% Nan Mask. (Not Fully Used Here.)
NanMask = false(size(sx));
sx(NanMask) = nan;
sy(NanMask) = nan;

%% Integration

tic
ztfli2 = tfli2(sx,sy,x,y);
toc

tic
zhfli2 = hfli2(sx,sy,x,y);
toc

tic
zsli2 = sli2(sx,sy,x,y);
toc

tic
zsli2p = sli2p(sx,sy,x,y);
toc


%% Show reconstruction errors.

eztfli2 = ztfli2 - z - mean(ztfli2(~NanMask)) + mean(z(~NanMask));
ezhfli2 = zhfli2 - z - mean(zhfli2(~NanMask)) + mean(z(~NanMask));
ezsli2  =  zsli2 - z - mean( zsli2(~NanMask)) + mean(z(~NanMask));
ezsli2p = zsli2p - z - mean(zsli2p(~NanMask)) + mean(z(~NanMask));

figure;
imshow([eztfli2, ezhfli2; ezsli2, ezsli2p],[]);
colormap jet;
colorbar;
