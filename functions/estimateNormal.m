function normal = estimateNormal(data, tree, query, radius, min_neighbors)
% Surface normal estimation using eigen decomposition

% Find neighbors within the fixed radius 
idx = rangesearch(tree, query, radius);                                                   
idxs = idx{1};
neighbors = [data(idxs(:),1) data(idxs(:),2) data(idxs(:),3)];
if size(neighbors) < min_neighbors
    normal = {1};
    return;
end
% Compute the covariance matrix C
C = cov(neighbors);
% Compute the eigenvector of C
[v, lambda] = eig(C);
% Find the eigenvector corresponding to the minimum eigenvalue in C
[~, i] = min(diag(lambda));
% Normalize
normal = v(:,i) ./ norm(v(:,i));
end 
