function iso_vectors = computeIsoVectors(vertices, edges, iso_points, vector_field)
    
    %% Version 3
    % This is used when the iso -points are computed using Marching
    % tetrahedra algorithm (i.e., iso-points must be on edges of tet mesh)
    
    % Note - compute_iso_vectors_2 is more generalised as it can be used
    % with an input point cloud set.

    % vertices: matrix (n, 3) containing the vertex coordinates
    % edges: matrix (m, 2) where each row gives the indices of the vertices that form an edge
    % iso_points: matrix (k, 3) containing the iso-point coordinates
    % vector_field: matrix (n, 3) containing the vector field values at each vertex
    % iso_vectors: matrix (k, 3) containing the interpolated vector at each iso-point
    
    num_iso_points = size(iso_points, 1); % Number of iso-points
    iso_vectors = zeros(num_iso_points, 3); % Initialize matrix to store the iso-vectors

    % Iterate through each iso-point
    for i = 1:num_iso_points
        iso_point = iso_points(i, :);
        found = false; % Flag to track if the edge containing the iso-point has been found
        
        % Search through all edges to find the edge containing the iso-point
        for e = 1:size(edges, 1)
            % Get the indices of the vertices that form the current edge
            v1_idx = edges(e, 1);
            v2_idx = edges(e, 2);
            
            % Get the coordinates of the two vertices of the edge
            V1 = vertices(v1_idx, :);
            V2 = vertices(v2_idx, :);
            
            % Compute the direction vector of the edge
            edge_dir = V2 - V1;
            
            % Compute the relative position t along the edge
            t = dot(iso_point - V1, edge_dir) / dot(edge_dir, edge_dir);
            
            % Check if the iso-point lies within the edge segment (0 <= t <= 1)
            if t >= 0 && t <= 1
                % Interpolate the vector field values at the edge vertices
                v1_vector = vector_field(v1_idx, :);
                v2_vector = vector_field(v2_idx, :);
                
                % Compute the interpolated vector at the iso-point
                iso_vectors(i, :) = (1 - t) * v1_vector + t * v2_vector;
                
                found = true;
                break; % Exit the loop since we found the containing edge
            end
        end
        
        % If no edge was found, display a warning
        if ~found
            warning('Iso-point %d does not lie on any edge.', i);
        end
    end
end
