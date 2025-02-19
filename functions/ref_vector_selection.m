function [best_candidates,min_variation_angles] = ref_vector_selection (orientation_vector_domains, local_max_indices,boundary_p)
    %Author : Joseph Jayakody
    %% Code - This is to find reference vectors at each local maxima point of the manufacturability map
    %% The input: Orientation vector domains (i.e., nozzle vector cones at each local maxima)
    %% For each candidate vector (i.e., each member of vector domain) at a local max point, we find
    %%      1. overhang using estimated boundary surface normal
    %%      2. variance w.r.t. to mean control vector
    %% Then, the best candidate vector is sorted as the one that generates minimum variance whilst adhereing
    %% to the overhang limit for support-free printing (i.e., 135 degrees). Repeated for every local max point.


    z_axis = [0  0  1]; %this can be involved with the rotated angle
    overhang_threshold = 140; variation_threshold = 60;
    min_variation_angle = 100;
    
    %% Note - variation threshold is to select the reference vectors from the nozzle cone
    % Overhang limit of 135 works with variation of 40
    % Overhang limit is 145 when variation is limited to 30

    %Estimated surface normal vectors at local max points
    est_sn_vectors = boundary_normal_vector_estimation (boundary_p, local_max_indices);     
    
    

    for f =1:size (orientation_vector_domains,1)

        current_sn_vector =  est_sn_vectors (f,:);

        current_domain = orientation_vector_domains {f};
        %zaxis_rep= repmat (z_axis, size(current_domain,1),1);

        mean_control_vector = mean (current_domain);

        for k=1:size(current_domain)
            current_candidate_vector = current_domain (k,:);
            overhang_angles(k,:)=atan2d(norm(cross(current_sn_vector,current_domain(k,:))), dot(current_sn_vector,current_domain(k,:)));
            
            % w.r.t. z-axis
            %variation_angles(k,:) = atan2d(norm(cross(z_axis,current_candidate_vector)), dot(z_axis,current_candidate_vector));
            %w.r.t. to mean control vector
            variation_angles(k,:) = atan2d(norm(cross(mean_control_vector,current_candidate_vector)), dot(mean_control_vector,current_candidate_vector));

            if  variation_angles(k,:) < min_variation_angle
                if variation_angles(k,:)<= variation_threshold && overhang_angles(k,:) <= overhang_threshold 
                    best_candidates {f,:} = current_candidate_vector;
                    min_variation_angles {f,:} = variation_angles(k,:);
                end
            end
        end
    end
end