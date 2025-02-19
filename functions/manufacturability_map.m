function [ local_maxima_values, local_maxima_indices, local_max_points, linear_map] = manufacturability_map (multifields,boundary_p,sigma)
    
    %% This map is a linear addition of the input constraint fields
    % Our goal is to find the local max points on the mesh geometry
    
    % Note - neighbourhood size parameter  is sigma (depends on mesh
    % quality)
    
    % User should define an appropriate size for the model - we use 10 for
    % most models, and 5 for others

    % Author - Joseph Jayakody
    linear_map = 0;
    for i=1:size(multifields,2)
        linear_map = linear_map + multifields (:,i);
    end

    
      
    neighbors = distancebased_vertex_neighborhood(boundary_p, sigma);
    
    % Find local maxima
    numVertices = size(boundary_p, 1);
    local_maxima_indices = [];
    
    for i = 1:numVertices
        is_local_max = true;
        
        current_value = linear_map(i);
        
        for j = 1:length(neighbors{i})
            neighbor_index = neighbors{i}(j);
            
            if linear_map(neighbor_index) > current_value
                is_local_max = false;
                break;
            end
            
        end
        if is_local_max
            local_maxima_indices = [local_maxima_indices, i];
        end
    end
    
    local_max_points = boundary_p ( local_maxima_indices,:);
    
    local_maxima_values = linear_map(local_maxima_indices);


end