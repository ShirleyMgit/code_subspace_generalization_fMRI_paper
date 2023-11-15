clear all 
close all
clc

mat_fig = [1 , 0.8, 0.3, 0.1; 0.8, 1, 0.1, 0.3; 0.3, 0.1, 1, 0.8; 0.1, 0.3, 0.8, 1];

figure(1)
imagesc(mat_fig)
xticks([])
yticks([])


voxels_illustration = randn(100,40);
figure(2)
imagesc(voxels_illustration)
xticks([])
yticks([])


