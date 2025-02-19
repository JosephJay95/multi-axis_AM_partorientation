function [Mesh_saliency_mutliscale, Mesh_saliency, local_max_points] = compute_saliency (boundary_p, F)
    
    %Author - Joseph Jayakody
    %
    
    % E should be 0.3% of the length of the bounding box  
    bbox_min = min(boundary_p);
    bbox_max = max(boundary_p);
    E = 0.03 * norm(bbox_max - bbox_min);

    k2 = discrete_mean_curvature(boundary_p,F);


    %% Saliency computation try - Chang Lee (2005)

    sigma_fine = 1*E; % Define fine scale sigma
    sigma_coarse = 2 * sigma_fine; % Coarse scale sigma
    
    % Step 2: 
    % Compute the neighborhoods for each scale
    % Here note that this function finds neighbourhood vertices based on a
    % threshold value (e.g., sigma_fine), similar to Piere's vring function
    
    vring_fine = distancebased_vertex_neighborhood(boundary_p, sigma_fine);
    vring_coarse = distancebased_vertex_neighborhood(boundary_p, sigma_coarse);
    
    % Step 3: Compute Gaussian weighted average mean curvature G[M,sigma]
    
    G_fine = gaussian_weighted_curvature(k2, boundary_p, vring_fine, sigma_fine);
    G_coarse = gaussian_weighted_curvature(k2, boundary_p, vring_coarse, sigma_coarse);
    
    % Step 5: Compute saliency
    Mesh_saliency = abs(G_fine - G_coarse);


    %% Step 6: Extend for multiple scales
    % You can plot saliency at different scales if you want!
    scales = E * (1:5);
    multi_scale_saliency = zeros(size(Mesh_saliency));
    
    for i = 1:length(scales)
        sigma_i = scales(i);
        vring_fine = distancebased_vertex_neighborhood(boundary_p, sigma_i);
        vring_coarse = distancebased_vertex_neighborhood(boundary_p, 2 *sigma_i);
        G_fine_i = gaussian_weighted_curvature(k2, boundary_p, vring_fine, sigma_i);
        G_coarse_i = gaussian_weighted_curvature(k2, boundary_p, vring_coarse, 2 * sigma_i);
        multi_scale_saliency = multi_scale_saliency + abs(G_fine_i - G_coarse_i);
    end
    
    Mesh_saliency_mutliscale = multi_scale_saliency / length(scales);
    
    %% Normalisation process
    % Set the scale parameter sigma for neighborhood selection
    sigma = 10;  
    neighbors = distancebased_vertex_neighborhood(boundary_p, sigma);
    
    % Find local maxima
    numVertices = size(boundary_p, 1);
    local_maxima_indices = [];
    
    for i = 1:numVertices
        is_local_max = true;
        
        current_value = Mesh_saliency_mutliscale(i);
        
        for j = 1:length(neighbors{i})
            neighbor_index = neighbors{i}(j);
            
            if Mesh_saliency_mutliscale(neighbor_index) > current_value
                is_local_max = false;
                break;
            end
            
        end
        if is_local_max
            local_maxima_indices = [local_maxima_indices, i];
        end
    end
    
    local_max_points = boundary_p ( local_maxima_indices,:);
    
    local_maxima_values = Mesh_saliency_mutliscale(local_maxima_indices);

    % Calculate global max M and average m of local maxima
    M = max(Mesh_saliency_mutliscale);
    m = mean(local_maxima_values);
    
    % Compute d and modify function values
    d = (M - m)^2;
    Mesh_saliency_mutliscale = Mesh_saliency_mutliscale * d;
    
    % Normalize the interpolated values to be within [0, 1]
    Mesh_saliency_mutliscale = Mesh_saliency_mutliscale / max(Mesh_saliency_mutliscale);


end