function visualize(this, t, options) 
% visualize cellLayer
%
% visualize(t)
% visualize(t, options)
%
% t:                time
%
% options:          
% -cellIndex:       show cell index
% -firstVertex:     color first vertex of the cell
% -transparent:     transparent faces, or 1 color transparent
% -edgeColor:       color of edges
% -lineWidth:       edge thickness
% -colorTable:      table of values for each cell to base coloring on
% -showBdryBonds:

nCells = numel(this.cells{t});

% define lattice variables lattmin style
cells = {this.cells{t}.bondInd};

maxNN = max(cellfun(@numel,cells));
cmp = lines(maxNN);
%maxNN = max(cellfun(@numel,cells)); 

bonds = cat(2,cat(1,this.bonds{t}.vertInd));%, cat(1,this.bonds{t}.cellInd));

if ~isempty(this.vertices2)
    verts = this.vertices2{t};
else
    verts = this.vertices{t};
end

% set options
if nargin == 2
    options = struct('cellIndex', 0, 'firstVertex', 0);
end

if isfield(options, 'edgeColor')
    edgeColor = options.edgeColor;
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
    
    if isfield(options, 'cellColorRange')
        
        minVal = options.cellColorRange(1);
        maxVal = options.cellColorRange(2);
    else
        
        minVal = min(CT);
        maxVal = max(CT);
    end
    
    cmap = jet(256);
    
    if isfield(options, 'colorbar') && options.colorbar
        colorbar('YTick',[0 255],...
         'YTickLabel',{num2str(minVal), num2str(maxVal)})
    end
    colorsIdx = round(255*min(CT - minVal, maxVal-minVal)/(maxVal-minVal) + 1);
    faceCol = cmap(colorsIdx,:);
    
% else contruct one from neighbor number
else
    faceCol = zeros([nCells 3]);
    
    for ci = 1 : nCells
        
        C = this.cells{t}(ci).bondInd;
        
        if ~isempty(C)
            
            vertsOfC = (bonds(C,1:2));
            vertsOfC = vertsOfC(:,1);
        
            % use black for cells with too many neighbors
            if length(vertsOfC) <= maxNN
                faceCol(ci,:) = cmp(length(vertsOfC),:);
            else
                faceCol(ci,:) = [0 0 0];
            end
        end
    end
end

% % visualize cells
% for ci = 1 : length(cells)
% 
%     if ~isempty(cells{ci})
%         
%         vertsOfC = (bonds(cells{ci}(:),1:2));
%         vertsOfC = vertsOfC(:,1);
%         
%         if isfield(options, 'transparent') && strcmp(options.transparent,'all')
%             
%             fcol = 'None';
%         else
%             if isfield(options, 'transparent') && isfield(options, 'colorTable')...
%                 && options.transparent == options.colorTable(ci)
%                 fcol = 'None';
%             else
%                 fcol = faceCol(ci,:);
%             end
%         end
%         patch(verts(vertsOfC,1), verts(vertsOfC,2), verts(vertsOfC,3), faceCol(ci,:),...
%                    'EdgeColor', edgeColor, 'FaceColor',fcol,'LineWidth',lw);
%     end
% end

faces = nan([nCells maxNN]);

for ci = 1:length(cells) 
    C = cells{ci};
    faces(ci,1:numel(C)) = bonds(C,1);
end

fv = struct();
fv.Vertices = verts;
fv.Faces = faces;

if isfield(options, 'transparent') && options.transparent
    fv.FaceVertexCData = 0*verts + 1;
else
    fv.FaceVertexCData = faceCol;
end

% if isfield(options, 'transparent') && options.transparent
%     edgeColorTable = verts*0; 
%     for c = 1:3
%         edgeColorTable(:,c) = edgeColorTable(:,c) + edgeColor(c);
%     end
%     patch('Faces',faces,'Vertices',verts,'FaceColor','none','EdgeColor','flat', 'FaceVertexCData', edgeColorTable)
% else
%     patch('Faces',faces,'Vertices',verts,'FaceVertexCData',faceCol)
% end

patch(fv)
%shading flat
shading faceted
axis equal

% visualize boundary bonds
if isfield(options, 'showBdryBonds') && options.showBdryBonds
    
    bdryBonds = this.getBdryBonds(t); % this will not work: placeholder
    
    for ci = 1:length(bdryBonds)

        %v1 = bonds(bdryBonds(ci),1);
        %v2 = bonds(bdryBonds(ci),2);

        v1 = bdryBonds(ci).vertInd(1);
        v2 = bdryBonds(ci).vertInd(2);

        hold on;
        quiver(verts(v1,1), verts(v1,2),...
            verts(v2,1)- verts(v1,1), verts(v2,2)-verts(v1,2), 0,...
                                    'Color', edgeColor, 'LineWidth', lw);
        hold off;
    end
end

% label cells by their index 
if isfield(options, 'cellIndex')
    
    nCells = numel(cells);
    if numel(options.cellIndex) == 1 && options.cellIndex
        cellIndex = 1:nCells;
    elseif numel(options.cellIndex) > 1
        cellIndex = options.cellIndex;
    else
        cellIndex = [];
    end
    
    if ~isempty(cellIndex)
        for ci = 1:nCells
            VofC = bonds(cells{ci}, 1:2);
            zpos = 0.01; % make sure text is above patch
            text(mean(verts(VofC,1)), mean(verts(VofC,2)), zpos,...
                                    num2str(cellIndex(ci)), 'Color', 'w');
        end
    end
end

% label first vertex
if isfield(options, 'firstVertex') && options.firstVertex
    hold on;
    for ci = 1:nCells
        VofC = bonds(cells{ci}, 1:2);
        quiver(verts(VofC(1),1), verts(VofC(1),2),...
            verts(VofC(2),1) - verts(VofC(1),1),...
            verts(VofC(2),2) - verts(VofC(1),2), 0, '.g', 'ShowArrowHead', 'on');
    end
    hold off;
end

end
