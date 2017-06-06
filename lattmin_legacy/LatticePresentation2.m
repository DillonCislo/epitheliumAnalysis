function LatticePresentation2(g, varargin) 
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
% -lineWidth:       edge thickness
% -colorTable:      table of values for each cell to base coloring on

maxNN = 14;
cmp = lines(maxNN);

if nargin == 2
    options = varargin{1};
else
    options = struct('cellIndex', 0, 'firstVertex', 0);
end

if isfield(options, 'edgeColor') && numel(options.edgeColor)==1
    
    edgeColor = options.edgeColor;
    
elseif isfield(options, 'edgeColor') && numel(options.edgeColor)>1
    
    edgeColor = [0 0 0];
    CT = options.edgeColor;
    
    if isfield(options, 'edgeColorRange')
        
        minVal = options.edgeColorRange(1);
        maxVal = options.edgeColorRange(2);
    else
        
        minVal = min(CT);
        maxVal = max(CT);
    end
   
	edgeCData = min(CT - minVal, maxVal-minVal)/(maxVal-minVal);
    
else
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
    minVal = -5;%-max(abs(CT))/2;%0;%min(CT);
    maxVal = 5;%max(abs(CT))/2;%max(CT)+1;
    CT(CT < minVal) = minVal;
    CT(CT > maxVal) = maxVal;
    cmap = jet(256);
    colorsIdx = round(255*min(CT - minVal, maxVal-minVal)/(maxVal-minVal) + 1);
    colors = cmap(colorsIdx,:);
    
% else contruct one from neighbor number
else
    colors = zeros([length(g.cells) 3]);
    
    for ci = 1 : length(g.cells)
        
        if ~isempty(g.cells{ci})
            
            verts = (g.bonds(g.cells{ci}(:),1:2));
            verts = verts(:,1);
        
            % use black for cells with too many neighbors
            if length(verts) <= maxNN
                colors(ci,:) = cmp(length(verts),:);
            else
                colors(ci,:) = [0 0 0];
            end
        end
    end
end

% visualize cells
% TODO: this would probably be faster using the fv structure
for ci = 1 : length(g.cells)

    if ~isempty(g.cells{ci})
        
        verts = (g.bonds(g.cells{ci}(:),1:2));
        verts = verts(:,1);

        if (isfield(options, 'transparent') && options.transparent) ||...
           (isfield(options, 'transparentColor') && options.transparentColor == CT(ci))
       
            patch(g.verts(verts,1), g.verts(verts,2), colors(ci,:),...
                    'EdgeColor',edgeColor,'FaceColor','None','LineWidth',lw);
        else
            patch(g.verts(verts,1),g.verts(verts,2), colors(ci,:),...
                'EdgeColor',edgeColor,'FaceColor', colors(ci,:),'LineWidth',lw);
        end
        
        if isfield(options,'edgeColor') && numel(options.edgeColor) > 1
            cdata = edgeCData(circshift(g.cells{ci},[0 0]),:);
            patch(g.verts(verts,1), g.verts(verts,2), cdata,...
                    'EdgeColor','flat','FaceColor','None','LineWidth',2,...
                    'CDataMapping','scaled');
        end
    end
end
caxis([0 1]);
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

% label cells by their index 
if isfield(options, 'cellIndex')
    
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
