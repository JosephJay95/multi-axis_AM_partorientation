function [orientation_vector_domains] = vector_domain_sampling (points)

    numSamples = 200;
    
    
    for k = 1:size (points)
        
        %% Computing vector domain for each point and its reference direction 
        % i.e., The vector domain is built around the reference direction.

        origin = points (k,:); % Replace with your 3D point coordinates
        direction = [ 0  0  1]; % Replace with your 3D vector components

        %direction = direction / norm(direction); % Normalize the direction vector
        
        %% Step 1: Sample points on the upper hemisphere
        Rx_samples = 2*pi*rand(numSamples, 1); % Azimuthal angle: 0 to 2*pi
        Ry_samples = pi/2*rand(numSamples, 1); % Polar angle: 0 to pi/2
            
        % Convert spherical coordinates to Cartesian coordinates
        x_samples = cos(Rx_samples) .* sin(Ry_samples);
        y_samples = sin(Rx_samples) .* sin(Ry_samples);
        z_samples = cos(Ry_samples);
        samples = [x_samples, y_samples, z_samples];       
        
        % Since the direction vector is [0, 0, 1], no need to rotate the samples

        % Step 2: Translate points to the origin
        translated_samples = samples + origin;
        
              
        
        %% Step 4: Rotate the reference vector along the directions of the sampled points
        
        % Compute the new reference vectors for each sample direction
        % Translate the sample directions to the origin to get the relative vectors
        sample_directions = translated_samples - origin;
        
        % Normalize the sample directions
        sample_directions = sample_directions ./ vecnorm(sample_directions, 2, 2);
                
        % For the fixed direction [0, 0, 1], the sample directions are already aligned
        % Simply use the sample directions as the rotated vectors

        orientation_vector_domains {k,:} = sample_directions;
    end

end
