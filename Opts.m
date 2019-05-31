% SuperVoxel multi-modal segmentation

% Introduction
% The algorithm used in this code is the implementation of our paper [1] 
% which was fulfilled in the Lab of Vision Engineering at the University of 
% Lincoln (UK). It is a modification of the method Simple Linear Iterative 
% Clustering (SLIC) which was proposed by Achanta et al. [2]. 
% Our method is optimized for medical images such as MRI, CT, etc. The 
% contributions of our codes compared to conventional 2D and 3D superpixel
% are as follows:
% •	Multi-modal input (works for single-modal, as well)
% •	Taking the spatial resolution of the medical images into account, i.e.
% the voxel resolution in X and Y directions and the slice thickness.


% parameters
Iterations = 10;
Compactness = 0.2;  % Compactness of the superpixel is equivalent to m
voxel_X = 10;
voxel_Y = 10;
voxel_Z = 4;

VoxelSize = [0.9375 0.9375 2.8];
term_x = 1/VoxelSize(1);
term_y = 1/VoxelSize(2);
term_z = 1/VoxelSize(3);