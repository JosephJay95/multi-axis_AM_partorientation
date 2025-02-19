function [final_score] = compliance_evaluation_optimisation (boundary_v_field, ref_vectors, boundary_p)
    % Author: Joseph Jayakody
    % Note - This is for optimisation. Because the objective just computes
    % the total penalty score (i.e., final score). So it will be minimised.

    % our penalty score is normalised;

    %% Note that our compliance evaluation (i.e., penalty score) is based on the 
    %% boundary tool orientation vector field. 

    z_axis = [0 0 1]; % This can be involved with the rotated angle
    overhang_threshold = 135;
    

    % Estimated boundary surface normal vectors at all boundary_p
    B_surface_normal = boundary_normal_vector_estimation_all(boundary_p);
    
    % Mean vectors of the smooth vector field and ref vectors - To find the variance
    Mean_vector1 = mean(boundary_v_field);
    Mean_vector2 = mean(ref_vectors);
    
    % Initialize arrays to store overhang and variation angles
    overhang_angles = zeros(size(boundary_v_field, 1), 1);
    variation_angles = zeros(size(boundary_v_field, 1), 1);
    
    % Initialize variables to accumulate penalties
    overhang_penalty = 0;
    variance_penalty = 0;

    % Calculate angles and penalties for each boundary point
    for m = 1:size(boundary_v_field, 1)
        % Overhang angle between smooth vectors and estimated surface normal vectors
        overhang_angle = atan2d(norm(cross(B_surface_normal(m,:), boundary_v_field(m,:))), dot(B_surface_normal(m,:), boundary_v_field(m,:)));
        overhang_angles(m) = overhang_angle;
        
        % Variation angle between smooth vectors and the z-axis
        variation_angle = atan2d(norm(cross(Mean_vector1, boundary_v_field(m,:))), dot(Mean_vector1, boundary_v_field(m,:)));
        variation_angles(m) = variation_angle;
        
        % Compute overhang penalty
        % >= might be important to make this differentiable
        if overhang_angle >= overhang_threshold
            overhang_penalty = overhang_penalty + (overhang_angle - overhang_threshold);
        end
        
        % Compute variance penalty
        variance_penalty = variance_penalty + variation_angle;
    end
    
    % USER DEFINED WEIGHTS
    W1 = 0.6;
    W2 = 0.4;

    % Total compliance score is the sum of overhang penalty and variance penalty
    penalty = (W1)*overhang_penalty + (W2)*variance_penalty;
    
    %normalising step
    penalty_score = penalty/size(boundary_p,1);
    
    %Compliance of the boundary tool vector field
    final_score = penalty_score;
end
