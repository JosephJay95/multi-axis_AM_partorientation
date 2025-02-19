function unique_edges = computeTetMeshEdges( T)
    % computeTetMeshEdges computes the unique edges of a tetrahedral mesh.
    % 
    % Inputs:
    %   V - Nx3 matrix of vertex coordinates
    %   T - Mx4 matrix of tetrahedral elements (each row contains 4 node indices)
    % 
    % Output:
    %   unique_edges - Px2 matrix of unique edges (each row contains 2 node indices)
    
    % Extract edges from each tetrahedron (6 edges per tetrahedron)
    edges = [T(:, [1 2]);  % Edge between node 1 and node 2
             T(:, [1 3]);  % Edge between node 1 and node 3
             T(:, [1 4]);  % Edge between node 1 and node 4
             T(:, [2 3]);  % Edge between node 2 and node 3
             T(:, [2 4]);  % Edge between node 2 and node 4
             T(:, [3 4])]; % Edge between node 3 and node 4

    % Sort the nodes in each edge to avoid duplicates
    sorted_edges = sort(edges, 2);

    % Remove duplicate edges
    unique_edges = unique(sorted_edges, 'rows');
end
