function transformedVertices = meshrotation_yaxis_only (vertices, sample)
        
        % Input - sample (between 0 to 1): 0 - null rotation, 1 - full rotation
        
        % This rotation is fixed along X or Y or both XY axes. So not
        % entirely unit sphere sampling However, we can extend this part at anytime..  
        % (refer to reference_vector selection code)
         
        % Author - Joseph Jayakody

        Ry = 2*pi*sample; % Range: 0 to 2*pi
        %Ry = pi*sample(2); % Range: 0 to pi
        centroid = mean(vertices);


        % Compute the rotation matrix around the x-axis (Rx)
%         Rx_matrix = [1, 0, 0, 0;
%                      0, cos(Rx), -sin(Rx), 0;
%                      0, sin(Rx), cos(Rx), 0;
%                      0, 0, 0, 1];
        
        % Compute the rotation matrix around the y-axis (Ry)
        Ry_matrix = [cos(Ry), 0, sin(Ry), 0;
                     0, 1, 0, 0;
                     -sin(Ry), 0, cos(Ry), 0;
                     0, 0, 0, 1];
        
        % Compute the transformation matrix for the current sample pair
        %T = Rx_matrix * Ry_matrix;
        T = Ry_matrix;
        
         % Apply the rotation to the mesh coordinates with respect to the original mesh  
        homogeneousVertices = [vertices - centroid, ones(size(vertices, 1), 1)];%adding 1s as 4th column
        transformedVertices = (T * homogeneousVertices')';
        
        % Adding back the centroid
        transformedVertices = transformedVertices(:, 1:3) + centroid;

end