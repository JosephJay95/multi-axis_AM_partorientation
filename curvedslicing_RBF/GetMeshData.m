function [Tri, X, fmt] = GetMeshData(TR)
% GetMeshData extracts faces and vertices from mesh in various formats

if isa(TR,'triangulation') || isa(TR,'TriRep')
    Tri = TR.ConnectivityList;
    X = TR.Points;
    fmt = 1;
elseif iscell(TR) && numel(TR)==2
    Tri = TR{1}; X = TR{2};
    fmt = 2;
elseif isstruct(TR) && isfield(TR,'faces') && isfield(TR,'vertices')
    Tri = TR.faces;
    X = TR.vertices;
    fmt = 3;
else
    error('Unrecognized mesh format');
end
end
