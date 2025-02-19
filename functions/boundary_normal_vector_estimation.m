function est_vectors = boundary_normal_vector_estimation (boundary_p, indices)

%% This code ONLY estimates boundary surface normal vectors at points defined by 'indices'.
%  boundary_p are the points used in 'esitmateNormal' function.

%Author - Joseph Jayakody

                BV= boundary_p;
                data = BV;
                tree = KDTreeSearcher(data);
                 
    
                % Below parameters (for surface normal estimation) can be changed
                radius = 5;
                min_neighbors = 4;
                          
                ref = mean (BV);

                for k= 1:size (BV)
                    query= BV(k,:);
                    est_surface_normal(k,:)=estimateNormal(data, tree, query, radius, min_neighbors);
                end

                
                B_surface_normal = est_surface_normal;
                vectors_to_surface = BV - ref;

                % Compute the dot product between each surface normal vector and the corresponding vector to the surface point
                dot_products = dot(B_surface_normal, vectors_to_surface, 2);

                % Find the indices of the inward-pointing surface normal vectors
                inward_indices = find(dot_products < 0);
                
                % Flip the orientation of the inward-pointing surface normal vectors
                B_surface_normal(inward_indices, :) = -1 * B_surface_normal(inward_indices, :);

                est_vectors = B_surface_normal (indices,:);
                est_vectors = est_vectors/ norm (est_vectors);

end