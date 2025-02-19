function B_surface_normal = boundary_normal_vector_estimation_all (boundary_p)

%% This code ALL estimates boundary surface normal vectors.
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

                

end