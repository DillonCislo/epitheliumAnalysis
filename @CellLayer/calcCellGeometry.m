function calcCellGeometry(this, t, varargin)
    % CALCCELLGEOMETRY Calculate geometric properties of cells
    %
    % calcCellGeometry(t)
    % calcCellGeometry(t, SOI, chartName, vectorPatch)
    %
    % SOI:          SurfaceOfInterest object
    % chartName:    chart in which cellLayer was made
    % vectorPatch:  reference vectors for orientation
    %
    % called without extra arguments for flat tissue geometry, 
    % with extra arguments for curved tissue
    
    %assert(~isempty(this.L), 'Label matrix CellLayer.L needs to be defined');
    
    if isempty(this.vertices2)
        V = this.vertices{t};
    else
        V = this.vertices2{t};
    end
    
    %----------------------------
    % flat case
	%----------------------------
    
	if nargin == 2
        
        disp('calculating flat cell geometries');

        % use regionprops to calculate some of them
        if ~isempty(this.L)
            
            stats = regionprops(this.L, 'Centroid', 'Area', 'MajorAxisLength',...
                                'MinorAxisLength', 'Orientation');
        else
            disp('implement version without label matrix');
        end
        
        for ci = 1:length(this.cells{t})

            l = this.cells{t}(ci).label;
            bondInd = this.cells{t}(ci).bondInd;
            
            if numel(bondInd) > 0
                naiveNeighbors = cat(1,this.bonds{t}(bondInd).cellInd);
                naiveNeighbors = naiveNeighbors(naiveNeighbors(:,2) > 0, 2);
                nn = naiveNeighbors;
            else
                nn = [];
            end

            geom = struct();
            geom.type = 'naive';
            
            % neighbor number
            geom.neighborNumber = numel(nn);
            
            % the region props ones
            if ~isempty(this.L)
                
                geom.centroid = stats(l).Centroid;
                geom.area = stats(l).Area;
                geom.majorAxisLength = stats(l).MajorAxisLength;
                geom.minorAxisLength = stats(l).MinorAxisLength;
                geom.orientation = stats(l).Orientation;
            end
%             
%             % circulation and tension
%             mech = this.cells{t}(ci).localMI();
%             geom.circulation = mech.circ;
%             geom.localT = mech.T(:,2);
%             geom.Tanisotropy = mech.Tanisotropy;
%             geom.Tprincipal = mech.Tprincipal;
%             geom.forceCov = mech.forceCov;
%             
            % set cell property
            this.cells{t}(ci).setGeometry(geom);

            % warning for less than 2 vertices 
            % for an epithelium this means some segmentation error, for a
            % culture with holes it can be ok
            if numel(bondInd) <= 2                    
                disp(['cell ' num2str(l) ' has 2 or less vertices!']);
            end
        end
    
    %----------------------------
	% curved case
    %----------------------------
    
    else
        
        disp('calculating metric corrected cell geometries');
                
        % assign input arguments
        SOI = varargin{1};
        chartName = varargin{2};
        vectorField = varargin{3};
        
%         % add 3D vertex positions
%         vertIdx = {this.vertices2{t}(:,1), this.vertices2{t}(:,2)};
%         verts3D = SOI.embedding.getPatch('cylinder1_index').apply(vertIdx);
%         this.vertices{t} = cat(2,verts3D{:});
%         
        % use regionprops to calculate flat properties, then correct using
        % metric
        stats = regionprops(this.L, 'Centroid', 'MajorAxisLength',...
                                'MinorAxisLength', 'Orientation');

        if isempty(SOI.g.patches)
            SOI.NCalcInducedMetric(t);
        end
                            
        % metric and metric determinant
        chart = SOI.atlas.getChart(chartName);
        gpatch = SOI.g.getPatch(chart.domain.name);
        metric = gpatch.getTransform(chartName).apply;
        
        detg = gpatch.determinant.apply{1};
        dA = sqrt(detg)*chart.image.stepSize(1)*chart.image.stepSize(2);
        
        assert(gpatch ~=0, 'calculate metric first, SOI.calcInducedMetric');

        for ci = 1:length(this.cells{t})
            
            l = this.cells{t}(ci).label;
            geom = struct();
            geom.type = 'proper';
            
            % neighbor number
            geom.neighborNumber = numel(this.cells{t}(ci).bondInd);
            
            % area 
            cidx = this.L == l;
            Acorrected = sum(dA(cidx)); %A = sum(cidx(:))
            geom.area = Acorrected;
            this.cells{t}(ci).setGeometry(geom);

            % warning for less than 2 vertices 
            % for an epithelium this means some segmentation error, for a
            % culture with holes it can be ok
            if geom.neighborNumber <= 2                    
                disp(['cell ' num2str(l) ' has 2 or less vertices!']);
            end
        end
    end
end
