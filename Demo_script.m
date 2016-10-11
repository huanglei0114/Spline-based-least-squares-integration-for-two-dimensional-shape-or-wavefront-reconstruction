% Demo_script.

%% Clear everyting existing.

clear ;
close all;
clc;

%% 0: Load Data.

load('Step_01_GradientData.mat','x','y','z','tsx','tsy','nsx','nsy');

% Choose the slope to integrate.
sx = nsx;
sy = nsy;

% NaN Mask. Not Fully Used Here.
NaNmask = false(size(sx));
sx(NaNmask) = NaN;
sy(NaNmask) = NaN;

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


%% Show the errors.

eztfli2 = ztfli2 - z - mean(ztfli2(~NaNmask)) + mean(z(~NaNmask));
ezhfli2 = zhfli2 - z - mean(zhfli2(~NaNmask)) + mean(z(~NaNmask));
ezsli2  =  zsli2 - z - mean( zsli2(~NaNmask)) + mean(z(~NaNmask));
ezsli2p = zsli2p - z - mean(zsli2p(~NaNmask)) + mean(z(~NaNmask));

figure;
imshow([eztfli2, ezhfli2; ezsli2, ezsli2p],[]);
