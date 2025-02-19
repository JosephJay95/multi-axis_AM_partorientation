    clear;close all;clc;
    addpath(genpath('Libraries'));  
      
    model = createpde(1);
    importGeometry(model,'new5_original.stl'); 
    b=generateMesh(model,'Hmax',2,'GeometricOrder','linear');
    V= b.Nodes';
    T=b.Elements';
    min_h=min(V(:,3));

    % Input - iso-points
    load ('knot_100_layers.mat');
    iso_points = cell2mat(ISO_PT2);
    iso_faces = cell2mat(ISO_CONN2);

    F = boundary_faces(T);
    boundary_p_indices=unique(F);
    boundary_p= V(boundary_p_indices,:);

    % Input - Saddle points and indices for Knot model at initial OR.
    %[saddle_p, saddle_idx] = compute_saddle_points (V,T);
    load('knot_saddle_points_and_indices.mat')
    saddle_points = saddle_p;
    saddle_id = saddle_idx;

    % Define the objective function - X axis only
    objective_func = @(x) part_orientation_optimiser_Xaxis_only(boundary_p, F, saddle_id, x);

    x0 = 0.2; % Example initial guess, we first tested with 0.5
    lb = 0; % Lower bound
    ub = 1; % Upper bound
    Iterations_no = 5;

    % Define the optimization options for pattern search
    options = optimoptions('patternsearch', 'Display', 'iter', 'PlotFcn', @psplotbestf, 'MaxIterations', Iterations_no );

    tic;
    [x_opt, fval] = patternsearch(objective_func, x0, [], [], [], [], lb, ub, options);    
    toc;
    

    optimal_orien_vertices = meshrotation_xaxis_only (boundary_p,x_opt);
    
    figure(1);
    scatter3(optimal_orien_vertices(:,1),optimal_orien_vertices(:,2),optimal_orien_vertices(:,3),'blue','filled');
    axis equal; axis off;



    %% ISO - VECTOR FIELD GENERATION
    
    % Scalar fields
    [TS_map_values, ~] = compute_topological_significance(optimal_orien_vertices,F, saddle_id);
    [Mesh_saliency_mutliscale, ~, ~] = compute_saliency (optimal_orien_vertices, F);
    [curvedness, ~, ~] = compute_local_shape (optimal_orien_vertices, F);
    
    % Salient feature map and points
    multifields = [curvedness  Mesh_saliency_mutliscale  TS_map_values];
    sigma = 10; 
    [ ~, local_maxima_indices, local_max_points] = manufacturability_map (multifields,optimal_orien_vertices, sigma);
    
    % Control vector selection
    [orientation_vector_domains] = vector_domain_sampling (local_max_points);
    [best_candidates,~] = ref_vector_selection (orientation_vector_domains, local_maxima_indices,optimal_orien_vertices);
    best_candidates = cell2mat(best_candidates);
    
    % Smooth boundary vector field design
    tic;
    A = cotmatrix (optimal_orien_vertices,F); 
    B = zeros (size(optimal_orien_vertices,1),3);      
    [smooth_v1 ,~ ]= min_quad_with_fixed (A,B,local_maxima_indices,best_candidates) ;
    best_candidates = best_candidates/ norm(best_candidates);
    toc;
    
    figure(1);
    quiver3(optimal_orien_vertices(:,1),optimal_orien_vertices(:,2),optimal_orien_vertices(:,3),smooth_v1(:,1),smooth_v1(:,2),smooth_v1(:,3), 'r','LineWidth',1, 'DisplayName', 'Boundary vector field'); hold on;
    quiver3(local_max_points(:,1),local_max_points(:,2),local_max_points(:,3),best_candidates(:,1),best_candidates(:,2),best_candidates(:,3), 'g','LineWidth',2, 'DisplayName', 'Best reference vectors'); 
    axis off; axis equal;

    % Volumetric vector field design
    A2 = cotmatrix (V ,T); 
    B2 = zeros (size(V ,1),3);        
    [smooth_v2 ,preA ]= min_quad_with_fixed (A2,B2,local_maxima_indices, best_candidates) ;
    
    % Iso-vector computation
    edges = computeTetMeshEdges(T);
    NEW_vectors = computeIsoVectors(V, edges, iso_points, smooth_v2);
