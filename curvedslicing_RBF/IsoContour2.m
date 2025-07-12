function [Contours, Levels, VData, Handles] = IsoContour(TR, F, Levels, showPlot)
% IsoContour: Extracts and plots isocontours (level sets) on a triangular mesh.
%
% INPUTS:
%   - TR        : Triangular mesh as triangulation, TriRep, or {faces, vertices}.
%   - F         : N x 1 scalar field at mesh vertices.
%   - Levels    : 1 x K vector of isovalues or scalar indicating # of levels.
%   - showPlot  : Logical or axis handle to show the result.
%
% OUTPUTS:
%   - Contours  : 1 x K cell array of Nx3 line segments per level.
%   - Levels    : 1 x K vector of isovalues used.
%   - VData     : 1 x K cell array of edge interpolation data.
%   - Handles   : 1 x K plot handles for contours.

% ————————————————————————
% 1. Input parsing
% ————————————————————————
if nargin < 2 || isempty(TR) || isempty(F)
    error('Must provide mesh and scalar field');
end

[Faces, Vertices, fmt] = GetMeshData(TR);
if fmt > 1
    TR = triangulation(Faces, Vertices);
end

F = F(:);
if numel(F) ~= size(Vertices,1)
    error('Length of F must match number of vertices');
end

% Handle isovalues
Fmin = min(F); Fmax = max(F);
if nargin < 3 || isempty(Levels)
    Levels = (Fmin + Fmax)/2;
elseif isscalar(Levels)
    if Levels == 1
        Levels = (Fmin + Fmax)/2;
    else
        dF = (Fmax - Fmin)/Levels;
        Levels = linspace(Fmin + dF/2, Fmax - dF/2, Levels);
    end
elseif numel(Levels) == 2 && Levels(1) == Levels(2)
    Levels = Levels(1);
end

if nargin < 4, showPlot = false; end
plotAxes = [];

% ————————————————————————
% 2. Initialize output
% ————————————————————————
nLevels = numel(Levels);
Contours = cell(1, nLevels);
VData = cell(1, nLevels);
Handles = nan(1, nLevels);

% ————————————————————————
% 3. Setup visualization
% ————————————————————————
if islogical(showPlot) && showPlot
    figure('Color','w');
    plotAxes = axes; hold on;
    trimesh(TR, 'Parent', plotAxes, ...
        'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8], 'FaceAlpha', 0.7);
    axis equal; lighting gouraud; camlight headlight;
elseif ishandle(showPlot)
    plotAxes = showPlot; hold(plotAxes, 'on');
end

% ————————————————————————
% 4. Process each level
% ————————————————————————
for i = 1:nLevels
    isoVal = Levels(i);
    segs = []; interpData = [];

    for triID = 1:size(Faces,1)
        verts = Faces(triID,:);
        fVals = F(verts);
        coords = Vertices(verts,:);

        % Find edges crossing the level
        below = fVals < isoVal;
        above = fVals > isoVal;

        if all(below) || all(above) || all(fVals == isoVal)
            continue; % no intersection or full triangle on contour
        end

        pairs = [1 2; 2 3; 3 1];
        crossings = [];

        for j = 1:3
            a = pairs(j,1);
            b = pairs(j,2);
            f1 = fVals(a); f2 = fVals(b);
            if (f1 - isoVal) * (f2 - isoVal) < 0 || f1 == isoVal || f2 == isoVal
                % Linear interpolation
                if f2 ~= f1
                    t = (isoVal - f1) / (f2 - f1);
                else
                    t = 0;
                end
                pt = (1 - t)*coords(a,:) + t*coords(b,:);
                crossings(end+1,:) = pt; %#ok<AGROW>
                interpData(end+1,:) = [verts(a), verts(b), t, 0]; %#ok<AGROW>
            end
        end

        % Add segment if two intersection points were found
        if size(crossings,1) == 2
            segs = [segs; crossings(1,:); crossings(2,:); NaN(1,3)]; %#ok<AGROW>
        end
    end

    Contours{i} = segs;
    VData{i} = interpData;

    % Plot
    if ~isempty(segs) && ~isempty(plotAxes)
        Handles(i) = plot3(plotAxes, segs(:,1), segs(:,2), segs(:,3), ...
            'k', 'LineWidth', 1.5);
    end
end

end
