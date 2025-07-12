% This code contains the RBF-weighted geodesic scalar field implementation
% for the curved layer slicing method proposed by the following paper [1]
% "Convexity and Surface Quality Enhanced Curved Slicing for Support-Free
% Multi-Axis Fabrication" - Jayakody et al. (2023)

% Notes: GA-based optimisation is not included. The user can adapt the
% current scalar field with their preferred optimisation model 

% Tool orientation planning can be implemented using the code related to
% the following paper [2] "A salient vector field-driven part orientation
% selection for multi-axis 3D printing" - Jayakody et al. (2025)

% Dependencies - GPtoolbox and functions related to [2] paper.


clear;close all;clc;
addpath(genpath('Libraries'));  
 
model = createpde(1);
importGeometry(model,'kitten_big.stl'); 
b=generateMesh(model,'Hmax',2,'GeometricOrder','linear');
V= b.Nodes';
T=b.Elements';

min_h=min(V(:,3));


tic;


F = boundary_faces(T);
boundary_p_indices=unique(F);
boundary_p= V(boundary_p_indices,:);

% Testing points produced from optimisation
load ('optimal_points_kitten.mat');
centroids = computeTetrahedronCentroids(T, V);

 
sigma = 10; % Spread parameter - this can be changed.

% % Define the  RBF function with Euclidean distances
gaussian_rbf = @(x, c, sigma) exp(-sum((x - c).^2, 2) / (2 * sigma^2));

% Initialize the interpolated values
interpolated_values_1 = zeros(size(T, 1), 1);

% Compute the RBF values and interpolate
for i = 1:size(all_points, 1)
    % Compute RBF values for each vertex in the mesh
    rbf_values = gaussian_rbf(centroids, all_points(i, :), sigma);
  
     % Accumulate the minimum RBF values
     interpolated_values_1 = interpolated_values_1 + rbf_values;
end

weights  = (1- ( interpolated_values_1  / max( interpolated_values_1 )));
weights = weights+ 0.05;


tol = min(V(:,3))+0.05;
[~, base_indices] = base_layer_extraction (V,tol);
[height_field_weighted,X]=heat_geodesics_original_weighted(V,T,base_indices, weights);

[height_field_weighted1,X1]=heat_geodesics_original(V,T,base_indices);


% Geodesic field
figure(40);
axis equal; grid off; axis off;
tsurf(F,V,'CData', height_field_weighted1); shading interp;
colorbar; axis equal;grid off; axis off;

TR={F,V};
ISO_H = (0:max(height_field_weighted1)/145:max(height_field_weighted1))'; 
[Q1,~,~,~]=IsoContour2(TR,height_field_weighted1,ISO_H,gca);

% Weighted field
figure(41);
axis equal; grid off; axis off;
tsurf(F,V,'CData', height_field_weighted); shading interp;
colorbar; axis equal;grid off;
 

TR={F,V};
ISO_H = (0:max(height_field_weighted)/145:max(height_field_weighted))'; 
[Q,~,~,~]=IsoContour2(TR,height_field_weighted,ISO_H,gca);
D = Q{77}; hold on;
scatter3 ( D(:,1),D(:,2),D(:,3),'red', 'filled'); 


 



