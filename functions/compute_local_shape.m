function [curvedness, SI, local_max_points] = compute_local_shape (boundary_p, F)
    %Author - Joseph Jayakody
    %
    K = discrete_gaussian_curvature(boundary_p,F);
    H = discrete_mean_curvature(boundary_p,F);
    M = massmatrix(boundary_p,F);
    k = H + [1 -1].*sqrt(H.^2 - M*K);
    
    % principal curvatures
    k1=k(:,1); k2=k(:,2);

    curvedness=real(sqrt((k1.^2+ k2.^2)/2))*10;

    SI=0;            
           
    for u=1:size(k1)
       if (k1(u)>k2(u))
           SI(u,:)=real((2/pi)*atan((k2(u)+k1(u))/(k2(u)-k1(u))));
       else
           SI(u,:)=0;
       end
    end

    %% Normalisation process
    % Set the scale parameter sigma for neighborhood selection
    sigma = 10;  
    neighbors = distancebased_vertex_neighborhood(boundary_p, sigma);
    
    % Find local maxima
    numVertices = size(boundary_p, 1);
    local_maxima_indices = [];
    
    for i = 1:numVertices
        is_local_max = true;
        
        current_value = curvedness(i);
        
        for j = 1:length(neighbors{i})
            neighbor_index = neighbors{i}(j);
            
            if curvedness(neighbor_index) > current_value
                is_local_max = false;
                break;
            end
            
        end
        if is_local_max
            local_maxima_indices = [local_maxima_indices, i];
        end
    end
    
    local_max_points = boundary_p ( local_maxima_indices,:);
    
    local_maxima_values = curvedness(local_maxima_indices);

    % Calculate global max M and average m of local maxima
    M = max(curvedness);
    m = mean(local_maxima_values);
    
    % Compute d and modify function values
    d = (M - m)^2;
    curvedness = curvedness * d;
    
    % Normalize the interpolated values to be within [0, 1]
    curvedness = curvedness / max(curvedness);

end