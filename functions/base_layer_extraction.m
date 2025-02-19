function [base_V, base_indices] = base_layer_extraction (V,tolerance)
    
    height= V(:,3);
    min_height=min(height);
    range = (height>=min_height) & (height<=tolerance); %outputs logic 1 or 0 for V within this range (a positive tolerance to extract faces)
    value=height(range);  %  find height values in the defined base plane range
    %min_value = min(value); max_value = max(value); %extremal values of the base plane range
    n=size(V,1); % no. of V 
    
    %base plane vertex indices finding
    
    for i =1:n
        x=V(i,:);                       % current vertex
        if ismember (x(:,3),value)     % checking if current vertex's 'z' value is in the base range
            base_V(i,:)=x;          % if it is, then this is a base vertex
            % finally we find the indice of the base (current) vertex for heat geodesic eq.
            base_indices(i,:) = find(ismember(V, base_V(i,:),'rows'));   
        end
    end
 
    %removing 0 values
    base_V( ~any(base_V,2), : ) = [];
    base_indices( ~any(base_indices,2), : ) = [];
end