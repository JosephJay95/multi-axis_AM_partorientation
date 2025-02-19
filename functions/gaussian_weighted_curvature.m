
function G = gaussian_weighted_curvature(k, V, vring, sigma)
    G = zeros(size(k));
    for i = 1:length(k)        
        neighbours = vring{i}; % neighbour vertex indices
        distances = sqrt(sum((V(neighbours, :) - V(i, :)).^2, 2));
        weights = exp(-distances.^2 / (2 * sigma^2));

        %Gaussian weighted average mean curvature
        G(i) = sum(weights .* k(neighbours)) / sum(weights);
    end
end