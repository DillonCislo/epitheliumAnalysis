function g = GLattConversion2(Cells, bulkVerts, varargin) 
    % Converts Voronoi-structures into GLattice structures
    %
    % g = GLattConversion2(Cells, bulkVerts)
    % g = GLattConversion2(Cells, bulkVerts, trim)
    % g = GLattConversion2(Cells, bulkVerts, bdryBonds, bdryVerts)
    %
    % bulkVerts:    nx3 array of bulk vertex coordinates
    %               bulk means part of a cell, not bdryBond
    % Cells:        {v1 .. vn} cell array of vertex indices
    % trim:         trim lattice to create boundary bonds, this removes
    %               outside cells with twofold vertices and only keeps the
    %               edges attached to other cells, useful for mech inverse
    % bdryBonds:    supply outside cells with bonds not belonging to any
    %               cell as boundary conditions (form below)
    % bdryVerts:    vertices going with the bdryBonds 
    % 
    %
    % g is a structure with fields:
    %
    %   verts:        nx3 array of vertex positions
    %   cells:        {b1 .. bn} cell array of bonds, bi index into bonds
    %                 bonds are ordered ccw, but matlab y<-> makes
    %                 that cw
    %   bonds:        nx4 array with rows [v1 v2 c1 c2]
    %                 vi, ci, respectively indices into verts and cells
    %                 bonds are directional, so pointing from v1 to
    %                 v2 with c1 to the left and c2 to the right
    %                 each bond is therefore part of only one cell,
    %                 with its antibond being part of the
    %                 neighboring cell
    % bdryBonds:    array with rows of the form [bi, ci, vi]
    %               bi is the index of a boundary bond into g.bonds, 
    %               ci the cell it is attached to, vi the vertex of cell
    % cti2ci        trimmed cell index to original cell idx
    %
    %---------------------------------------------------------------------
    % Copyright 2014 Idse Heemskerk and Sebastian Streichan
    %---------------------------------------------------------------------
    
    if nargin ==4 
        
        bdryBonds = varargin{1};
        bdryVerts = varargin{2};
        
    elseif nargin == 3
        
        trim = varargin{1};
        bdryBonds = [];
        bdryVerts = [];
        
    elseif nargin == 2
        
        warning('GLattConversion2: trim lattice by default!');
        
        trim = 1;
        bdryBonds = [];
        bdryVerts = [];
    end

    NCells = length(Cells);
    NBulkVerts = length(bulkVerts);
    NBdryVerts = length(bdryVerts);
    g = struct('cells',[],'bonds',[],'verts',[]);
    
    % trivial indexing without trimming, overwritten for trimming case
    g.cti2ci = 1:NCells;
        
    % Populate Vertices
    g(1).verts = zeros(NBulkVerts + NBdryVerts, 3);
    g(1).verts(1:NBulkVerts, :) = bulkVerts; % was 1:2, WHY?
    if ~isempty(bdryVerts)
        g(1).verts(NBulkVerts + 1: end, :) = bdryVerts; % was 1:2
    end
    
    % Make sure cells are sorted! 
    % btw: clockwise is anti-clockwise in images because y is down
    % matlab sucks
    for ci = 1:NCells
        if ~isempty(Cells{ci})
            
            [~,Index] = CCWSortCell(Cells{ci},bulkVerts);
            c = Cells{ci};
            Cells{ci} = c(Index);
            
            % some code later I guess assumes this extra copy of the verst
            % vertex at the end of the cell
            Cells{ci}(length(Index)+1) = c(Index(1));
        end
    end
    
    % Convert cell entries into bonds and note: 
    %------------------------------------------------
    % convention for bond=[v1 v2 c1 c2], looking from v1 to v2 -> c1= left cell
    % and c2->right cell
    [g(1).cells, g(1).bonds] = BondsConversion(Cells,NBulkVerts);        
    
    %-------------------------------------------
    % boundary bonds from argument
    %-------------------------------------------
    
    % include boundary bonds
    if ~isempty(bdryVerts)
        
        NBulkBonds = size(g(1).bonds, 1);
        NBdryBonds = size(bdryBonds, 1);

        % indexing of bdry vertices is shifted in g.verts wrt bdryVerts
        bdryBonds(:,2) = bdryBonds(:,2) + NBulkVerts;
        % bdry bonds have no cell (0) on either side
        bdryBonds = [bdryBonds, zeros(size(bdryBonds))];
        g(1).bonds = cat(1, g(1).bonds, bdryBonds);

        % bdryBonds = [b c v]
        g(1).bdryBonds = zeros([size(bdryBonds,1) 3]);
        g(1).bdryBonds(:,1) = NBulkBonds+1:NBulkBonds+NBdryBonds;

        % now go through bonds in the cell and identify cell and vertex that a
        % boundary bond attaches to
        for ci = 1:NCells

            cedges = g.cells{ci};

            for bi = 1:length(cedges)

                % if two consecutive bonds have border the outside
                % the vertex in between connects a boundary edge
                if g.bonds(g.cells{ci}(bi), 4) == 0 &&...
                    g.bonds(g.cells{ci}(mod(bi, length(cedges))+1), 4) == 0

                    vout = g.bonds(g.cells{ci}(bi), 2);
                    bdryBondIdx = find(g.bonds(g.bdryBonds(:,1),1) == vout);
                    g.bdryBonds(bdryBondIdx, 2) = ci;
                    g.bdryBonds(bdryBondIdx, 3) = vout;
                end
            end
        end
        
    %-------------------------------------------
    % else trim lattice to make boundary bonds 
    %-------------------------------------------
    
    elseif isempty(bdryVerts) && trim

        % trimmed list of cells and index into it
        ctrimmed = {};
        cti = 1;
        cutbonds = [];
        cti2ci = [];

        for ci = 1:NCells
            
            % neigboring cell indices of bonds
            nci = g.bonds(g.cells{ci},4);
            
            % if the cell has two consecutive bonds on the outside, it will
            % have a twofold vertex and should be removed
            % bonds of the removed cells are stored for the next step
            if any(nci==0 & circshift(nci,-1)==0)
                cutbonds = [cutbonds, g.cells{ci}];
            else
                ctrimmed{cti} = g.cells{ci};
                cti2ci(cti) = ci;
                cti = cti + 1;
            end
        end

        g.cti2ci = cti2ci;

        % reverse indexing: if cell in trimmed bulk, gives index otherwise
        % zero
        ci2cti = zeros([NCells 1]);
        ci2cti(cti2ci) = 1:length(cti2ci);

        % bonds remaining in the trimmed bulk
        bulkBondIdx = unique([ctrimmed{:}]);
        
        % bonds cell indices need to be updated to trimmed cell indices
        for bi = 1:size(g.bonds,1)
            if g.bonds(bi,3)~=0
                g.bonds(bi,3) = ci2cti(g.bonds(bi,3));
            end
            if g.bonds(bi,4)~=0
                g.bonds(bi,4) = ci2cti(g.bonds(bi,4));
            end
        end

        % vertices in the bulk
        bulkVertIdx = unique(g.bonds(bulkBondIdx,1:2));        
        % quick lookup is vertex in bulk
        vertInBulk = false([size(g.verts,1) 1]);
        vertInBulk(bulkVertIdx) = 1;

        % bonds attached to cells in the trimmed list but not part of them
        % are now boundary bonds
        rmBonds = false([size(g.bonds,1) 1]);
        bdryBondsTmp = [0 0 0];
        bbi = 1;

        for bi = 1:length(cutbonds)
                
            % if a cut bonds has no vertex in the bulk (is not attached to
            % a remaining cell) we can remove it
            % if both vertices are in the bulk, it is a bulk antibond and
            % can also be removed
            if  ~any(vertInBulk(g.bonds(cutbonds(bi),1:2))) ||...
                all(vertInBulk(g.bonds(cutbonds(bi),1:2)))
                
                rmBonds(cutbonds(bi)) = 1;
            else
                
                vert1 = g.bonds(cutbonds(bi),1);
                vert2 = g.bonds(cutbonds(bi),2);

                % for outside bonds, I want the bulk vertex first
                if vertInBulk(vert1) 
                    
                    % determine bulk cell containing the bulk vertex
                    bulkBondNeighbor = g.bonds(g.bonds(:,1) == vert1 &...
                                                    g.bonds(:,3) ~= 0,:);
                    ci = bulkBondNeighbor(3);   
                    bdryBondsTmp(bbi, :) = [cutbonds(bi) ci vert1];
                    bbi = bbi + 1;
                    
                else    
                    if ~any(g.bonds(cutbonds,1) == vert2 &...
                            g.bonds(cutbonds,2) == vert1)
                        
                        %flip bond
                        g.bonds(cutbonds(bi),2) = vert1;
                        g.bonds(cutbonds(bi),1) = vert2;
                        
                        %determine bulk cell containing the bulk vertex
                        bulkBondNeighbor = g.bonds(g.bonds(:,1) == vert2 &...
                                                        g.bonds(:,3) ~= 0,:);
                        ci = bulkBondNeighbor(3);   
                        bdryBondsTmp(bbi, :) = [cutbonds(bi) ci vert2];
                        bbi = bbi + 1;
                    else
                        
                       rmBonds(cutbonds(bi)) = 1;
                    end
                end
            end

%                 hold on;
%                 quiver(g.verts(vi,1), g.verts(vi,2),...
%                     g.verts(ve,1)-g.verts(vi,1), g.verts(ve,2) -g.verts(vi,2),0,'.g');
%                 scatter(g.verts(vi,1), g.verts(vi,2),'.g');
%                 hold off;
        end

        % now some garbage collection
        % first remove bonds that were not needed and reindex them
        NremBonds = size(g.bonds,1) - sum(rmBonds);
        newBondIdx(~rmBonds) = 1:NremBonds;
        for ci = 1:numel(ctrimmed)
            ctrimmed{ci} = newBondIdx(ctrimmed{ci});
        end
        for bi = 1:size(bdryBondsTmp,1)
            bdryBondsTmp(bi,1) = newBondIdx(bdryBondsTmp(bi,1));
        end
        g.bonds(rmBonds, :) = [];
        g.cells = ctrimmed;
        g.bdryBonds = bdryBondsTmp;

        % then remove unused vertices and reindex
        vertsUsed = unique(g.bonds(:,1:2));
        vertsUsedIdx = zeros([size(g.verts,1) 1]);
        vertsUsedIdx(vertsUsed) = 1:length(vertsUsed);
        g.verts(vertsUsedIdx==0,:) = [];
        for bi = 1:size(g.bonds,1)
            g.bonds(bi,1:2) = vertsUsedIdx(g.bonds(bi,1:2));
        end
        for bi = 1:size(g.bdryBonds,1)
            g.bdryBonds(bi,3) = vertsUsedIdx(g.bdryBonds(bi,3));
        end
       
    else
        
        g.bdryBonds = [];
    end

    %-------------------------------------------
    % bond conversion
    %-------------------------------------------
    
    function [cells,bonds] = BondsConversion(VorCells,NVertices)
        
        % determine NBonds 
        NBonds = 0;
        for i = 1 : NCells
            NBonds = NBonds + length(VorCells{i}) - 1;
        end
 
        %error('bla');
        % the same 2 vertices can be connected in the same direction by
        % different cells if a 2-sided cell is involved
        % we store this in duplicates
        S=sparse(NVertices,NVertices);
        dups = [];
        for i = 1 : NCells
            
            Nsides = length(VorCells{i})-1; %(last vertex is duplicate)
            for j = 1 : Nsides             
                
                v1 = VorCells{i}(j);
                v2 = VorCells{i}(j+1);

                if Nsides > 2
                    S(v1,v2) = i;
                else
                    %warning(['dumping 2 sided cell ' num2str(i)]);
                    dups = cat(1, dups, [v1 v2 i]);
                end
            end
        end
        % NOTE: if there is more than one dup, 3 cell share a bond, this
        % can happen if a non-convex cell has the vertices order in the
        % wrong way
        
        bonds = zeros(NBonds,4);
        Counter = 1;
        cells = cell(1,NCells);
        
        for i = 1 : NCells
            
            cells{i} = zeros(1,length(VorCells{i})-1);

            for j = 1 : (length(VorCells{i})-1)
                
                v1 = VorCells{i}(j);
                v2 = VorCells{i}(j+1);
                
                bonds(Counter,:) = full([v1 v2 i S(v2,v1)]);
 
                cells{i}(j) = Counter;
                Counter = Counter + 1;
            end
        end

        % now deal with dups
        % these come from two sided cells, which happen on the image edge
        % or colony edge, there are then two bonds connecting the same
        % vertices in the same direction, the one belonging to the two
        % sided cell should have outside cell 0, while the other one should
        % have the 2-sided cell as its outside
        for i = 1:size(dups,1)
            
            cellidx = dups(i,3);
            
            % bond and antibond index
            bj = bonds(:,1) == dups(i,1) & bonds(:,2) == dups(i,2);
            abj = bonds(:,2) == dups(i,1) & bonds(:,1) == dups(i,2);
            
            neighbor = setdiff(bonds(bj, 3), cellidx);

            if numel(neighbor) > 1
                error('what?');
            end
            
            if numel(neighbor)==1 % change to ~isempty(neighbor) to find errors related to 4-folds etc
                
                bonds(bj & bonds(:,3)==neighbor, 4) = cellidx;
                bonds(abj & bonds(:,3)==cellidx, 4) = neighbor;
            end
        end
    end
end