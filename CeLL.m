classdef CeLL < handle
    % representing a cell in a CellLayer
    %
    % CeLL, because Matlab has defined Cell
    
    %---------------------------------------------------------------------
    % properties
    %---------------------------------------------------------------------    

    properties (Access = private)
        
        cellLayer
        time 
    end
    
    properties (SetAccess = protected)
        
        label
        
        bondInd
        
        state
        geometry        % structure containing geometric properties
        geometryProper
    end
    
    properties (Dependent = true)
        
        outside
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------    

    methods
         
        %------------------------------------------------------
        % constructor
        %------------------------------------------------------
        
        function this = CeLL(cellLayer, time, label, bondInd, varargin)
            % initialize CeLL
            % 
            % CeLL(cellLayer, time, label, bondInd)
            % CeLL(cellLayer, time, label, bondInd, stateNames)
            %
            % cellLayer:    CellLayer to which the CeLL belongs
            % bondInd:      List of bond indices
            % label:        Index in array of cellLayer.

            
            % to make arrays of an object, it needs to be possible to
            % call the constructor without arguments. matlab sucks
            if nargin==0
                return;
            end
            
            this.cellLayer = cellLayer;
            this.bondInd = bondInd;
            this.time = time;
            this.label = label;
        end 
        
        %------------------------------------------------------
        % get neighbors
        %------------------------------------------------------
        
        function [nn, nnInd] = getNeighbors(this)
            % GETNEIGHBORS return a cell array of neighboring cell objects
            %
            % [nn nnInd] = getNeighbors()
            %
            % nn:       array of CeLL objects
            % nnInd:    indices of neighbor cells in ccw order
            % 
            % this does not include the outside, use CeLL.outside to check
            
            bonds = this.cellLayer.getBond(this.time, this.bondInd);
            nnInd = cat(1,bonds.cellInd);
            nnInd = nnInd(:,2);
            
            % exclude outside from neighbors
            nn = this.cellLayer.getCell(this.time, nnInd(nnInd~=0));
        end
        
        function outside = get.outside(this)
            % ISOUTSIDE determine if cell is on outside of cell layer
            
            % bond index is not bond label!
            bonds = this.cellLayer.bonds{this.time}(this.bondInd);
            nn = [bonds.cellInd];
            if any(nn == 0)
                outside = true;
            else
                outside = false;
            end
        end
        
        %------------------------------------------------------
        % get verts
        %------------------------------------------------------
        
        function cv = getVerts(this)
            % GETNEIGHBORS return an array of the cell's vertex labels
            % Will be ordered.
            
            bonds = this.cellLayer.getBond(this.time, this.bondInd);
            cellVerts = [bonds.vertInd];
            cv = zeros(1,length(cellVerts)/2);
            v = cellVerts(1:2);
            cv(1) = v(ismember(v,cellVerts(end-1:end)));
            for i = 2:length(cv)
                v = cellVerts(2*i-1:2*i);
                cv(i) = v(ismember(v,cellVerts(2*i-3:2*i-2)));
            end

        end
        
        %------------------------------------------------------
        % setters
        %------------------------------------------------------
        
        function setBondInd(this, bondInd)
            % needed for relabeling after cleaning up bonds
            % I prefer it to making new CeLL objects
            
            %warning('use setBondInd with care');
            this.bondInd = bondInd;
        end
            
            
            
        function setState(this, stateVector)
            % SETSTATE set the state vector
            %
            % setState(stateVector)
            %
            % stateVector length should match stateNames
            
            if isempty(this.cellLayer.cellState)
                error('state variables were never defined');
                
            elseif numel(stateVector) ~= numel(this.cellLayer.cellState)
                
                error('state vector length does not match names length');
            else
                
                this.state = stateVector;
            end
        end
        
        function setLabel(this, label)
            % SETSTATE set the label (into CellLayer.L)
            %
            % setLabel(label)
            %
            % label:    index into label matrix
            
            this.label = label;
        end
        
        function setGeometry(this, geometry)
            % SETGEOMETRY set structure with geometric properties
            % 
            % setGeometry(geometry)
            %
            % these properties are more efficiently computed at the level
            % of CellLayer and then set by this method

            if strcmp(geometry.type,'naive')
                this.geometry = geometry;
            elseif strcmp(geometry.type,'proper')
                this.geometryProper = geometry;
            else
                warning('type unknown, did nothing');
            end
        end
        
        function n = getNormal(this)
            % GETNORMAL get cell normal
            %
            % n = getNormal()
            
            bonds = this.bondInd;
            
            t = this.time;
            CL = this.cellLayer;
            
            vertIdx = cat(1, CL.bonds{this.time}(bonds).vertInd);
            vertIdx = vertIdx(:,1);
            
            vertsOfCell = CL.vertices{t}(vertIdx,:);
            vertsOfCell2 = circshift(vertsOfCell, 1);

            bondVecs = vertsOfCell2 - vertsOfCell;
            circBonds = circshift(bondVecs(:,:),1);

            normal = cross(bondVecs, circBonds);
            for ci = 1:size(normal,1)
                normal(ci,:) = normal(ci,:)/norm(normal(ci,:));
            end
            
            n = mean(normal);
        end
        
        function mech = localMI(this)
            % LOCALMI local mechanical inverse
            % 
            % mech = localMI()
            %
            % mech:     structure

            l = this.label;
            bonds = this.bondInd;
            nn = numel(bonds);
            
            t = this.time;
            CL = this.cellLayer;
            V = CL.vertices2{t};
            
            vertIdx = cat(1, CL.bonds{this.time}(bonds).vertInd);
            vertIdx = vertIdx(:,1);
            
            % circulation
            circ = 0;
            
            % tensions for left, right and outgoing edges of vertex
            T = zeros([nn 3]);
            
            % unit vectors around the cell
            rijhat = zeros([nn 2]);
            
            % for each vertex
            for vi = 1:nn

                imin = mod(vi-2,nn) + 1;
                iplus = mod(vi,nn) + 1;

                leftVidx = vertIdx(imin);
                rightVidx = vertIdx(iplus);
                outVidx = setdiff(CL.getVertNeighbors(t, vertIdx(vi)), vertIdx(:));

                if numel(outVidx) > 1

                    disp(['higher fold vertex in cell ' num2str(l)]);
                    circ = NaN;
                    
                elseif isempty(outVidx)
                    
                    disp(['2-fold vertex in cell label ' num2str(l)]);
                else

                    % three outgoing unit vectors
                    rleft = V(vertIdx(vi),:) - V(leftVidx,:);
                    rright = V(vertIdx(vi),:) - V(rightVidx,:);
                    rout = V(vertIdx(vi),:) - V(outVidx,:);
                    
                    rleft = rleft/norm(rleft);
                    rright = rright/norm(rright);
                    rout = rout/norm(rout);

                    % add to rijhat for anisotropy calculation later
                    rijhat(vi,:)= rright(1:2);
                    
                    % solve for right and out tensions in terms of left
                    % we can write force balance as r_{i\alpha}T_\alpha = 0
                    % where i is space and alpha is edges out of the vertex
                    % so r_{ia}T_a = -r_{i1}T_1 with a=2,3
                    % so T_a = -r_{ia}^{-1} r_{i1} T_1 
                    % let R = r_ia
                    T(vi,:) = [1, (-inv([rright(:,1:2)' rout(:,1:2)'])*rleft(:,1:2)')'];
                    
                    % actually I want unit tension on the first cw
                    if vi==1 
                    	T(1, :) = T(1,:)./T(1,2);
                    end

                    % rescale to get constency
                    if vi > 1
                        T(vi,:) = T(vi,:)*T(vi-1,2);
                    end
                    
                    % calculate contribution to circulation
                    ratioOfSines = -cross(rleft, rout)/cross(rright, rout);
                    circ = circ + log(ratioOfSines);
                    
                    if ratioOfSines < 0
                        disp(['imaginary part in cell label ' num2str(l)]);
                    end
                end
            end
            
            %-----------------------------------------------------
            % tension anisotropy per cell
            %-----------------------------------------------------

            if numel(this.bondInd) <= 2
                disp(['cell ' num2str(this.label) ' has 2 or less neighbors!']);
            else
                % ANISOTROPY FROM FORCE COVARIANCE
                
                % tensions on the edges
                Tij = T(:,2);

                if any(abs(Tij) == Inf) || any(isnan(Tij))
                    disp(['cell ' num2str(ci) ' has Inf or Nan tension on some edge']);
                    Tanisotropy = NaN;
                else
                    % forces
                    Fij = rijhat.*[Tij Tij];
                    forceCov = cov(Fij);
                    [V,D] = eig(forceCov);

                    evals = max(D);
                    [majoreig, idx] = max(evals);
                    minoreig = min(evals);
                    majoraxis = V(:,idx);

                    if all(imag(evals) == 0)
                        Tanisotropy = sqrt(majoreig / minoreig) - 1;
                        cellPrincipalT = majoraxis';
                    else
                        disp('FIXME: imaginary or infinite anisotropy');
                        Tanisotropy = NaN;
                    end
                end
            end
            
            mech.T = T;
            mech.circ=  circ;
            mech.forceCov = forceCov;
            mech.Tanisotropy = Tanisotropy;
            mech.Tprincipal = cellPrincipalT;
        end
    end
end