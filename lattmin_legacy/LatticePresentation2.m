function LatticePresentation2(g, options) 
% display Lattmin lattice 
%
% LatticePresentation2(g)
% LatticePresentation2(g, options)
%
% g:                Lattmin structure
% options:          
% -cellIndex:       show cell index
% -firstVertex:     color first vertex of the cell
% -transparent:     transparent faces
% -edgeColor:       color of edges
% -edgeColorRange:  color range of lattice edges
% -lineWidth:       edge thickness
% -colorTable:      table of values for each cell to base coloring on
% -cellColorRange:  color range of lattice faces
% -includeEdgeColorbar: include a colorbar for edge based fields
% -includeFaceColorbar: include a colorbar for face based fields
% -edgeColormap:    colormap to use for edge based fields
% -cellColormap:    colormap to use for face based fields

maxNN = 14;
cmp = lines(maxNN);
cmp = [cmp; 0 0 0];

if (nargin < 2), options = struct('cellIndex', 0, 'firstVertex', 0); end

if isfield(options, 'edgeColor') && numel(options.edgeColor)==1
    
    % Single number index into a saved colormap
    % This option does not appear to be used...
    edgeColor = options.edgeColor;
    
elseif isfield(options, 'edgeColor') && numel(options.edgeColor)>1
    
    % Scalar data per edge with optional color re-scaling
    edgeColor = [0 0 0];
    CT = options.edgeColor;
    
    if isfield(options, 'edgeColorRange')
        
        minVal = options.edgeColorRange(1);
        maxVal = options.edgeColorRange(2);

    else
        
        minVal = min(CT);
        maxVal = max(CT);

    end

    % Average edge data for bond-antibond pairs (only one color can be
    % displayed per edge)
    bonds = sort(g.bonds(:, 1:2), 2);
    dupl = find_duplicate_rows(bonds);

    duplIDx = struct2cell(dupl);
    duplIDx = duplIDx(2,:).';
    assert(all(cellfun(@numel, duplIDx) == 2), 'Non-manifold edge found');
    duplIDx = [duplIDx{:}].';

    antiBondIDx = repmat((1:size(g.bonds,1)).', 1, 2);
    antiBondIDx(duplIDx, 1) = [duplIDx(:,2); duplIDx(:,1)];

    CT = (CT + CT(antiBondIDx(:,2))) ./ 2;

    rmIDx = antiBondIDx(:,1) > antiBondIDx(:,2);
    bonds(rmIDx, :) = [];
    CT(rmIDx, :) = [];

	% edgeCData = min(CT - minVal, maxVal-minVal)/(maxVal-minVal);

    CT(CT < minVal) = minVal;
    CT(CT > maxVal) = maxVal;

    if isfield(options, 'edgeColormap')
        edge_cmap = options.edgeColormap;
    else
        edge_cmap = parula(256);
    end

    valsN = round(((CT - minVal) ./ (maxVal-minVal)) .* 255)+1;
    valsN(isnan(valsN)) = 1;

    edgeColors = edge_cmap(valsN,:);
    edge_crange = [minVal maxVal];
    
else

    % Default to black edges
    edgeColor = [0 0 0];

end

if isfield(options, 'lineWidth')
    lw = options.lineWidth;
else
    lw = 1;
end

% if there is a color table, use it
if isfield(options, 'colorTable')

    CT = options.colorTable(:);

    if isfield(options, 'cellColorRange')

        minVal = options.cellColorRange(1);
        maxVal = options.cellColorRange(2);

    else
        
        minVal = min(CT);
        maxVal = max(CT);
    end

    CT(CT < minVal) = minVal;
    CT(CT > maxVal) = maxVal;

    if isfield(options, 'cellColormap')
        cell_cmap = options.cellColormap;
    else
        cell_cmap = parula(256);
    end

    valsN = round(((CT - minVal) ./ (maxVal-minVal)) .* ...
        (size(cell_cmap,1)-1))+1;
    valsN(isnan(valsN)) = 1;

    cellColors = cell_cmap(valsN,:);
    cell_crange = [minVal maxVal];
    
% else contruct one from neighbor number
else

    % NOTE: This assumes no boundary antibonds are maintained
    bdyBondIDx = find((g.bonds(:,3) == 0) | (g.bonds(:,4) == 0));

    cellColors = cellfun(@(x) numel(x) - sum(ismember(x, bdyBondIDx)), ...
        g.cells).';
    cellColors(cellColors > maxNN) = maxNN + 1;
    cellColors = cmp(cellColors, :);
    cell_cmap = cmp;
    cell_crange = [1 maxNN];

    % cellColors = zeros([length(g.cells) 3]);
    % 
    % for ci = 1:length(g.cells)
    %     
    %     if ~isempty(g.cells{ci})
    %         
    %         verts = (g.bonds(g.cells{ci}(:),1:2));
    %         verts = verts(:,1);
    %     
    %         % use black for cells with too many neighbors
    %         if length(verts) <= maxNN
    %             cellColors(ci,:) = cmp(length(verts),:);
    %         else
    %             cellColors(ci,:) = [0 0 0];
    %         end
    %     end
    %   
    % end

end

if ~isfield(options, 'transparent')
    options.transparent = false;
end

% Transparent face choice overrides face color bar option
if (~isfield(options, 'includeFaceColorbar') || ...
        options.transparent)
    options.includeFaceColorbar = false;
end

% Face color bar option overrides edge color bar option
if (~isfield(options, 'includeEdgeColorbar') || ...
        options.includeFaceColorbar)
    options.includeEdgeColorbar = false;
end


% Convert the cell connectivity list to a (NaN-padded)
% matrix Useful for plotting purposes
maxFaceSize = max(cellfun(@numel, g.cells));
cellFace = nan(numel(g.cells), maxFaceSize);
% bondv = g.bonds(:, 1:2);
for i = 1:size(cellFace,1)
    cellFace(i, 1:numel(g.cells{i})) = g.bonds(g.cells{i}, 1).';
end

% Visualize cells
dispCells = ~cellfun(@isempty, g.cells).';
if options.transparent

    patch('Faces', cellFace(dispCells, :), 'Vertices', g.verts(:, 1:2), ...
        'FaceColor', 'none', 'EdgeColor', edgeColor, 'LineWidth', lw);

else

    patch('Faces', cellFace(dispCells, :), 'Vertices', g.verts(:, 1:2), ...
        'FaceVertexCData', cellColors, 'FaceColor', 'flat', ...
        'EdgeColor', edgeColor, 'LineWidth', lw);

end

if (isfield(options, 'edgeColor') && numel(options.edgeColor) > 1)

    hold on
    patch( 'Faces', reshape([1:size(bonds,1),1:size(bonds,1)*2],[],3), ...
        'Vertices', [g.verts(bonds(:), 1), g.verts(bonds(:), 2)], ...
        'FaceVertexCData', repmat(edgeColors, 2, 1), ...
        'EdgeColor', 'interp' );
    hold off

end

if options.includeFaceColorbar
    cb = colorbar;
    colormap(cell_cmap);
    set(gca, 'Clim', cell_crange);
end

if options.includeEdgeColorbar
    cb = colorbar;
    colormap(edge_cmap);
    set(gca, 'Clim', edge_crange);
end


% OLD SLOW VERSION WITH MULTIPLE PATCH OBJECTS ****************************
% TODO: this would probably be faster using the fv structure
% for ci = 1 : length(g.cells)
% 
%     if ~isempty(g.cells{ci})
%         
%         verts = g.bonds(g.cells{ci}(:), 1:2);
%         verts = verts(:,1);
% 
%         if (isfield(options, 'transparent') && options.transparent) ||...
%            (isfield(options, 'transparentColor') && options.transparentColor == CT(ci))
%        
%             patch(g.verts(verts,1), g.verts(verts,2), cellColors(ci,:),...
%                     'EdgeColor',edgeColor,'FaceColor','None','LineWidth',lw);
%         else
%             patch(g.verts(verts,1),g.verts(verts,2), cellColors(ci,:),...
%                 'EdgeColor',edgeColor,'FaceColor', cellColors(ci,:),'LineWidth',lw);
%         end
%         
%         if isfield(options,'edgeColor') && numel(options.edgeColor) > 1
%             cdata = edgeCData(circshift(g.cells{ci},[0 0]),:);
%             patch(g.verts(verts,1), g.verts(verts,2), cdata,...
%                     'EdgeColor','flat','FaceColor','None','LineWidth',2,...
%                     'CDataMapping','scaled');
%         end
%     end
% end
% caxis([0 1]);
% *************************************************************************

axis equal

% visualize boundary bonds
for ci = 1:length(g.bdryBonds)
    
    v1 = g.bonds(g.bdryBonds(ci),1);
    v2 = g.bonds(g.bdryBonds(ci),2);
    
    hold on;
    quiver(g.verts(v1,1), g.verts(v1,2),...
        g.verts(v2,1)- g.verts(v1,1), g.verts(v2,2)-g.verts(v1,2), 0,...
                                'Color', edgeColor, 'LineWidth', lw);
    hold off;
    
end

if ~isfield(options, 'cellIndex')
    options.cellIndex = false;
end

% label cells by their index 
if options.cellIndex
    
    nCells = numel(g.cells);
    if numel(options.cellIndex) == 1 && options.cellIndex
        cellIndex = 1:nCells;
    elseif numel(options.cellIndex) > 1
        cellIndex = options.cellIndex;
    else
        cellIndex = [];
    end
    
    if ~isempty(cellIndex)
        for ci = 1:nCells
            VofC = g.bonds(g.cells{ci}, 1:2);
            zpos = 0.01; % make sure text is above patch
            text(mean(g.verts(VofC,1)), mean(g.verts(VofC,2)), zpos,...
                                    num2str(cellIndex(ci)), 'Color', 'k');
        end
    end
end

% label first vertex
if isfield(options, 'firstVertex') && options.firstVertex
    hold on;
    nCells = numel(g.cells);
    for ci = 1:nCells
        VofC = g.bonds(g.cells{ci}, 1:2);
        quiver(g.verts(VofC(1),1), g.verts(VofC(1),2),...
            g.verts(VofC(2),1) - g.verts(VofC(1),1),...
            g.verts(VofC(2),2) - g.verts(VofC(1),2), 0, '.g', 'ShowArrowHead', 'on');
    end
    hold off;
end

set(gcf,'Color','white');

end

function [ dupl, C ] = find_duplicate_rows( A )
%FIND_DUPLICATE_ROWS Finds duplicate rows in a matrix

[C, ia, ~] = unique( A, 'rows' );

if size(A,1) == size(C,1)
    % disp('There are no duplicate rows!');
    dupl = {};
    return;
end

rep_idx = setdiff(1:size(A,1), ia);
rep_val = unique( A(rep_idx,:), 'rows');

dupl_val = cell( size(rep_val,1), 1 );
dupl_idx = cell( size(rep_val,1), 1 );

for i = 1:size(rep_val,1)
    
    dupl_val{i} = rep_val(i,:);
    
    diffA = A - repmat( rep_val(i,:), size(A,1), 1 );
    dupl_idx{i} = find( sqrt(sum(diffA .* conj(diffA), 2)) < eps );
    
end

dupl = struct( 'val', dupl_val, 'idx', dupl_idx );

end

