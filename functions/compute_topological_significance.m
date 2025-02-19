function [output, local_max_points, local_maxima_indices] = compute_topological_significance(boundary_p,faces, reeb_idx)
    %Author - Joseph Jayakody%
    %We use Geodesics for for Gaussian RBF computation on the mesh (Praun et al. 2000 -vector field design)
    
    % Parameters for Gaussian RBF
    sigma = 7.5;

    % Initialize the interpolated values
    TS_map_values = zeros(size(boundary_p, 1), 1);
    
    % Compute the RBF values and interpolate
    for i = 1:size(reeb_idx, 1)
        
        % Compute geodesic distances from the current imp_point to all vertices
        [D, ~] = heat_geodesics_singlesource(boundary_p, faces, reeb_idx(i));
    
        % Compute RBF values using geodesic distances
        rbf_values = exp(-D.^2 / (2 * sigma^2));

        % Accumulate the interpolated values
        TS_map_values = TS_map_values + rbf_values;
    end
    
    output = TS_map_values;

    
    %% Normalisation process
    % Set the scale parameter sigma for neighborhood selection
    sigma = 10;  
    neighbors = distancebased_vertex_neighborhood(boundary_p, sigma);
    
    % Find local maxima
    numVertices = size(boundary_p, 1);
    local_maxima_indices = [];
    
    for i = 1:numVertices
        is_local_max = true;
        
        current_value = TS_map_values(i);
        
        for j = 1:length(neighbors{i})
            neighbor_index = neighbors{i}(j);
            
            if TS_map_values(neighbor_index) > current_value
                is_local_max = false;
                break;
            end
            
        end
        if is_local_max
            local_maxima_indices = [local_maxima_indices, i];
        end
    end
    
    local_max_points = boundary_p ( local_maxima_indices,:);
    
    local_maxima_values = TS_map_values(local_maxima_indices);

    % Calculate global max M and average m of local maxima
    M = max(TS_map_values);
    m = mean(local_maxima_values);
    
    % Compute d and modify function values
    d = (M - m)^2;
    TS_map_values = TS_map_values * d;
    
    % Normalize the interpolated values to be within [0, 1]
    TS_map_values = TS_map_values / max(TS_map_values);

end