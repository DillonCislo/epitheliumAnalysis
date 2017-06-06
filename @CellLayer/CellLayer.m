classdef CellLayer < handle_light
    % Represents a layer of cells as it evolves over time 
    % Keeps track of topology, geometry and dynamics
    
    %---------------------------------------------------------------------
    % properties
    %---------------------------------------------------------------------
    
    properties (SetAccess = protected) %(Access = private) 
        
        cells           % cell array {time}, CeLL object arrays
        bonds           % cell array {time}, Bond object arrays
        vertices2       % cell array {time}, N x 3 arrays of xy vertex coordinates
        vertices        % cell array {time}, N x 3 arrays of xyz vertex coordinates
        nTimePts        % total number of frames
    end

    properties (SetAccess = protected)
        
        cellState       % cell state variable names, cell array of strings
        bondState       % bond state variable names, cell array of strings
        L               % label matrix for bonds and cells
    end
    
    properties (Dependent = true)
        
        bdryBonds       % the outward bonds providing BCs, not part of cell
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------    

    methods
         
        %------------------------------------------------------
        % constructor and initialization
        %------------------------------------------------------
        
        function this = CellLayer(nTimePts, varargin)
            % initialize cell layer
            %
            % CellLayer(nTimePts)
            % CellLayer(nTimePts, cellState, bondState)
            %
            % nTimePts:     total number of frames
            % cellState:    cell array of strings labeling state variables
            %
            % initialize data using initTime or init

            this.cells = {};
            this.bonds = {};
            this.vertices2 = {};
            this.vertices = {};
            this.nTimePts = nTimePts;
            
            this.L = []; % label matrix has 3rd dimension time, 
                         % more efficient than cell array
            
            if nargin == 3
                this.cellState = varargin{1};
                this.bondState = varargin{2};
            elseif nargin == 1
                this.cellState = {};
            else
                error('incorrect number of arguments');
            end
        end 
        
        %------------------------------------------------------
        % cell, bond, vertex getters
        %------------------------------------------------------
        
        function c = getCell(this, t, cellLabel)
            % GETCELL return cell at time
            %
            % c = getCell(t,cellLabel)
            % 
            % t:            time
            % cellLabel:    cell label or list of cell labels
            % c:            CeLL object
            
            assert(nargin == 3, 'provide time and cell labels');
            
            cellLabels = [this.cells{t}.label];
            cellIdx = false(size(cellLabels));
            for i = 1:numel(cellLabel)
                cellIdx = cellIdx | (cellLabels == cellLabel(i));
            end

            c = this.cells{t}(cellIdx);
        end
        
        function cellNN = getCellNeighbors(this, t, cellLabel)
            % GETCELLNEIGHBORS return nearest neighbors of cell
            %
            % cellNN = getCellNeighbors(t, cellLabel)
            %
            % t:            time
            % cellLabel:    cell label
            % cellNN:       labels of nearest neighbors
            
            cellObj = this.getCell(t, cellLabel);
            
            if isempty(cellObj)
                cellNN = NaN;
            else
                bIdx = cellObj.bondInd;
                bCells = cat(1,this.bonds{t}(bIdx).cellInd);
                cellNNIdx = bCells(:,2);
                cellNNIdx = cellNNIdx(cellNNIdx > 0);
                cellNN = [this.cells{t}(cellNNIdx).label];
            end
        end
        
        function b = getBond(this, t, bondLabel)
            % GETBOND return bond at time
            %
            % b = getBond(t,i)
            % 
            % t:            time
            % i:            Bond label or list of bond labels
            % b:            Bond object
            
            assert(length(this.bonds) >= t && ~isempty(this.bonds{t}),...
                    'bonds at this time have not been initialized');
            
            bondLabels = cat(1,this.bonds{t}.label);
            bondIdx = false(size(bondLabels));
            for i = 1:numel(bondLabel)
                bondIdx = bondIdx | (bondLabels == bondLabel(i));
            end

            b = this.bonds{t}(bondIdx);
        end
        
        function vertNN = getVertNeighbors(this, t, vertLabel)
            % GETVERTNEIGHBORS returns all nearest neighbors of vertexLabel
            %
            % vertNN = getVertNeighbors(t, v)
            %
            % t:            time
            % v:            Vertex label
            % vertNN:       array of vertex labels
            
            bondVerts = cat(1,this.bonds{t}.vertInd);
            vertNN = bondVerts( bondVerts(:,1) == vertLabel |...
                                bondVerts(:,2) == vertLabel, :);
            vertNN = setdiff(unique(vertNN(:)), vertLabel);
        end
        
        function vertCellInd = getVertCells(this, t, vertLabel)
            % GETVERTCELLS return the cells of which a vertex is part
            %
            % vertCellInd = getVertCells(t, vertLabel)
            % 
            % t:            time
            % vertLabel:    Vertex label
            % cellLabels:   neighboring cell labels
            
            % find the bonds containing the vertex 
            bondVertInd = cat(1,this.bonds{t}.vertInd);
            vertBonds = bondVertInd(:,1) == vertLabel |...
                        bondVertInd(:,2) == vertLabel;

            bondCellInd = cat(1,this.bonds{t}.cellInd);
            vertCellInd = unique(bondCellInd(vertBonds,:));
            vertCellInd(vertCellInd==0)=[];
        end
        
        function bdryBonds = getBdryBonds(this, t)
            % GETBDRYBONDS return outward pointing bonds on the boundary
            %
            % bdryBonds = getBdryBonds(t)
            %
            % t:            time
            % bdryBonds:    array of outward pointing bonds providing BCs
            
            bondCellInd = cat(1, this.bonds{t}.cellInd);
            bdryBonds = this.bonds{t}(bondCellInd(:,1) == 0 & bondCellInd(:,2) == 0);
        end
        
        %------------------------------------------------------
        % old lattice getters
        %------------------------------------------------------
        
        function g = getLattmin(this, t)
            % GETLATTMIN return a Lattmin structure
            %
            % g = getLattmin(t)
            %
            % t:    time
            % g:    Lattmin structure
            
            % cells
            g = struct();
            g.cells = {this.cells{t}.bondInd};
            
            % bonds
            bondVertInd = cat(1, this.bonds{t}.vertInd);
            bondCellInd = cat(1, this.bonds{t}.cellInd);
            g.bonds = [bondVertInd, bondCellInd];
            
            % verts
            warning('distinction between 2D and 3D vertices is sloppy');
            if ~isempty(this.vertices2)
                g.verts = this.vertices2{t};
            else
                g.verts = this.vertices{t};
            end
            
            % boundary bonds : [bi vi ci]
            bbInd = bondCellInd(:,1) == 0 & bondCellInd(:,2) == 0;
            bbVertInd = bondVertInd(bbInd, 1);
            
            bbCellInd = zeros(size(bbVertInd));
            for i = 1:size(bbCellInd,1)    
                vertCells = this.getVertCells(t, bbVertInd(i));
                bbCellInd(i) = vertCells(1);
                if length(vertCells) > 1
                    disp('warning: bdrybond on 4-fold vertex');
                end
            end
            bbInd = find(bbInd);
            g.bdryBonds = [bbInd, bbCellInd, bbVertInd];
        end
        
        %------------------------------------------------------
        % other getters
        %------------------------------------------------------
        
        function array = getCellGeometry(this, t, property)
            % GETCELLGEOMETRY Return a cell geometry property as array
            %
            % T = getCellGeometry(t, property)
            % 
            % t:        time
            % property: string with property name
            % T:        array with indexing matching cell indexing 
            
            assert(isnumeric(t), 'first argument should be time index');
            cellGeometries = cat(1, this.cells{t}(:).geometry);
            array = cat(3, cellGeometries.(property));
            array = squeeze(array);
        end
        
        function array = getBondGeometry(this, t, property)
            % GETBONDGEOMETRY Return a bond geometry property as array
            %
            % T = getBondGeometry(t, property)
            %
            % t:        time
            % property: string with property name
            % T:        array with indexing matching bond indexing 
            
            bondGeometries = cat(1, this.bonds{t}(:).geometry);
            array = cat(3, bondGeometries.(property));
            array = squeeze(array);
        end
        
        function bondIdx = getStateBoundary(this, t, state, value)
            
            si = find(strcmp(this.cellState(:),state));
            bondCells = cat(1, this.bonds{t}.cellInd);
            
            nBonds = numel(this.bonds{t});
            RcellState = false([nBonds 1]);
            LcellState = false([nBonds 1]);
            insideB = bondCells(:,1) ~= 0 & bondCells(:,2) ~= 0;
            
            RcellState(insideB) = cat(1,this.cells{t}(bondCells(insideB,1)).state);
            LcellState(insideB) = cat(1,this.cells{t}(bondCells(insideB,2)).state);
            
            bondIdx = [ RcellState >= value & LcellState < value,...
                        RcellState < value & LcellState >= value ];
        end
        
        %------------------------------------------------------
        % setters
        %------------------------------------------------------

        function setCellLabels(this, t, labelVector)
            % SETCELLLABELS sets the label property of cells pointing into
            % a label matrix
            %
            % setCellLabels(t, labelVector)
            %
            % t:            time
            % labelVector:  vector containing cell labels
            
            for i = 1:length(this.cells{t})
                this.cells{t}(i).setLabel(labelVector(i));
            end
        end

        function setCellState(this, t, statrix)
            % SETCELLSTATE sets the state property of cells 
            %
            % setCellLabels(t, stateVector)
            %
            % t:            time
            % statrix:      matrix of cell states
            %               Ncells x Nstates
            
            for i = 1:length(this.cells{t})
                this.cells{t}(i).setState(statrix(i,:));
            end
        end

        function setLabelMatrix(this, t, L)
            % SETLABELMATRIX set label matrix to which cell labels refer
            %
            % setLabelMatrix(t, L)
            %
            % t:    time
            % L:    label matrix
            
            % if L didn't exist, make it now as a single matrix for all
            % times (for more efficient memory allocation)
            if isempty(this.L)
                this.L = zeros(size(L,1), size(L,2), this.nTimePts, 'uint16');
            end
            
            this.L(:,:,t) = L;
        end
        
        %------------------------------------------------------
        % remove cells
        %------------------------------------------------------
        
        function removeCells(this,t,cellIdx)
            % removeCells(t,cellIdx)
            
            %error('bla');
            this.cells{t}(cellIdx) = [];
            bondsUsed = unique(cat(2,this.cells{t}(:).bondInd));
            
            nBonds = numel(this.bonds{t});
            newNBonds = numel(bondsUsed);
            relabel = zeros([nBonds 1]);
            relabel(bondsUsed) = 1:newNBonds;
            
            this.bonds{t} = this.bonds{t}(bondsUsed);
            nCells = numel(this.cells{t});
            for i = 1:nCells
                this.cells{t}(i).setBondInd(relabel(this.cells{t}(i).bondInd)');
            end
        end
        
        %------------------------------------------------------
        % state dependent label matrix
        %------------------------------------------------------
        
        function specialL = cellL(this, t, index)
            % return the label matrix for cells only
            
            cellLabels = cat(1,this.cells{t}.label);
            
            if nargin == 3
                cellLabels = cellLabels(index);
            end
            
            labelIntervals = idx2interval(cellLabels);
            nIntervals = size(labelIntervals,1);
            
            mask = false(size(this.L(:,:,t)));
            
            for i = 1:nIntervals
                mask = mask | ( this.L(:,:,t) >= labelIntervals(i,1) &...
                                    this.L(:,:,t) <= labelIntervals(i,2) );
            end
            
            specialL = this.L(:,:,t);
            specialL(~mask) = 0;
        end

        function specialL = bondL(this, t, index)
            % return the label matrix for bonds only
            
            bondLabels = cat(1,this.bonds{t}.label);
            bondLabels(bondLabels == 0) = [];
            
            if nargin == 3
                bondLabels = bondLabels(index);
            end

            labelIntervals = idx2interval(bondLabels);
            nIntervals = size(labelIntervals,1);
            
            mask = false(size(this.L(:,:,t)));
            
            for i = 1:nIntervals
                mask = mask | ( this.L(:,:,t) >= labelIntervals(i,1) &...
                                    this.L(:,:,t) <= labelIntervals(i,2) );
            end
            
            specialL = this.L(:,:,t);
            specialL(~mask) = 0;
        end
        
        %------------------------------------------------------
        % dual lattice
        %------------------------------------------------------
        
        function dualCL = dualize(this)
            % construct the dual lattice
            
            % looks funny in visualize because vertices are not ordered
            % cellLayer2 = dualLayer.dualize;
            % cellLayer2.visualize;
            % this also doesn't work

            dualCL = CellLayer(this.nTimePts);
            
            for t = 1:this.nTimePts
                
                nCells = numel(this.cells{t});
                verts = this.vertices{t};
                bondVert = cat(1,this.bonds{t}.vertInd);
                bondCell = cat(1,this.bonds{t}.cellInd);

                % dual vertices are cells with cell center of mass as position
                dualVerts = zeros([nCells 3]);
                for i = 1:nCells
                    dualVerts(i,:) = mean(verts(bondVert(this.cells{1}(i).bondInd,1),:));
                end

                % edges are dual to edges
                % in this case interior edges
                interiorBonds = bondCell(:,1)~=0 & bondCell(:,2)~=0;
                dualBondVert = bondCell(interiorBonds,:);
                dualBondCell = bondVert(interiorBonds,[2 1]);
                
                % relabel dual cells
                dualCellInd = unique(dualBondCell(:,1));
                nDualCells = numel(dualCellInd);
                relabel = zeros([max(dualCellInd) 1]);
                relabel(dualCellInd) = 1:nDualCells;
                dualBondCell = relabel(dualBondCell);

                % construct dual cells
                dualCells = {};
                n = 1;
                for i = 1:nDualCells
                    indualcell = find(dualBondCell(:,1) == i)';
                    % only if it is a triangle or more
                    if numel(indualcell) > 2
                        dualCells{n} = indualcell;
                        n = n+1;
                    end
                end
                nDualCells = n - 1;
                
                % order the bonds
                for  ci = 1:nDualCells
                    X = cat(1, dualBondVert(dualCells{ci},:));
                    P = [1 0 0];
                    for i = 2:numel(dualCells{ci})
                        P(i) = find(X(:,1) == X(P(i-1),2));
                    end
                    dualCells{ci} = dualCells{ci}(P);
                end

                dualg = struct();
                dualg.verts = dualVerts;
                dualg.bonds = [dualBondVert dualBondCell];
                dualg.cells = dualCells;
                dualg.bdryBonds = [];
                
                dualg
                
                dualCL.initTime(t, 'Lattmin', dualg);
            end
        end
    end    
end