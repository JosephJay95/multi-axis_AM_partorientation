function vring = distancebased_vertex_neighborhood(V, sigma)
    % This function finds neighbouring vertices based on the threshold value - sigma
    % Therefore the neighbourhood size is controlled by sigma
    % The operation is similar to Piere's vring function - returns neighbourhood as a cell
    %% Author - Joseph Jayakody 
    
    num_vertices = size(V, 1);
    vring = cell(num_vertices, 1);
    distances = pdist2(V, V); % Compute pairwise distances between all vertices

    for i = 1:num_vertices
        vring{i} = find(distances(i, :) < sigma);
    end
end