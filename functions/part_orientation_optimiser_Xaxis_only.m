 function [cost] = part_orientation_optimiser_Xaxis_only (boundary_p,F, reeb_idx, x)
%function [best_candidates, local_maxima_indices, cost] = part_orientation_optimiser_Xaxis_only (boundary_p,F, reeb_idx, x)

%  Note- that for fmincon algorithms, the function must only output a scalar value

    %% Objective function for the optimisation model - X axis rotation only
    % Author - Joseph Jayakody
 
    %% Input: x is a vector that stores rotation angles along X axis
    %% x is between 0 to 1, 0 - null rotation, 1- full rotation.
     
    global bestCandidates bestLocalMaxPoints 

    %x(1) = Rx, x(2) = Ry
    V_transformed = meshrotation_xaxis_only (boundary_p, x);
    
    % Manufacturability map
    [TS_map_values, ~] = compute_topological_significance(V_transformed,F, reeb_idx);
    [Mesh_saliency_mutliscale, ~, ~] = compute_saliency (V_transformed, F);
    [curvedness, ~, ~] = compute_local_shape (V_transformed, F);

    % Initial manufacturability map
    multifields = [curvedness  Mesh_saliency_mutliscale  TS_map_values];
    sigma=10;
    [ ~, local_maxima_indices, local_max_points] = manufacturability_map (multifields,V_transformed,sigma);
    
    [orientation_vector_domains] = vector_domain_sampling (local_max_points);
    [best_candidates,~] = ref_vector_selection (orientation_vector_domains, local_maxima_indices,V_transformed);
    best_candidates = cell2mat(best_candidates);

    %% Smooth vector field design
    A = cotmatrix (V_transformed,F); % contagent laplacian of the surface layer
    B = zeros (size(V_transformed,1),3);
    
    % we can actually use the pre computation in our optimisation loop
    [smooth_v1 ,~ ]= min_quad_with_fixed (A,B,local_maxima_indices,best_candidates) ;
    best_candidates = best_candidates/ norm(best_candidates);
    

    penalty_score  = compliance_evaluation_optimisation (smooth_v1, best_candidates, V_transformed); 
   
    % our optimisation model minimises the penalty cost (to improve compliance)
    % so the best part orientation has highest compliance (i.e., lowest cost) 

    cost  = penalty_score; 

 
end