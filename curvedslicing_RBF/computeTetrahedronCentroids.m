function centroids = computeTetrahedronCentroids(T, V)
    
    %% Author- Joseph Jayakody 
    % T: nx4 matrix - tetrahedra with each row containing the vertex indices
    % V: mx3 matrix - representing the vertex coordinates

    % Initialize an empty array to store the centroids
    centroids = zeros(size(T, 1), 3);

    % Loop through each tetrahedron
    for i = 1:size(T, 1)
        % Extract vertex indices of the current tetrahedron
        vertexIndices = T(i, :);

        % Extract coordinates of the vertices using the indices and reshape
        vertices = reshape(V(vertexIndices, :), [4, 3]);

        % Compute the centroid using the mean function along each dimension
        centroid = mean(vertices);

        % Store the centroid in the array
        centroids(i, :) = centroid;
    end
end