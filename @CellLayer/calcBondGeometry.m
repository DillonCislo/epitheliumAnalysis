function calcBondGeometry(this, t, varargin)
    % CALCBONDGEOMETRY Calculate geometric properties of cells
    %
    % calcBondGeometry(t)
    % calcBondGeometry(t, SOI, chartName, vectorPatch)
    %
    % SOI:          SurfaceOfInterest object
    % chartName:    chart in which cellLayer was made
    % vectorPatch:  reference vectors for orientation
    %
    % called without extra arguments for flat tissue geometry, 
    % with extra arguments for curved tissue
    
    %assert(~isempty(this.L), 'Label matrix CellLayer.L needs to be defined');
 
    %----------------------------
    % flat case
	%----------------------------
    
	if nargin == 2
        
        disp('calculating bond geometries');

        % use regionprops to calculate some of them
        if ~isempty(this.L)
            
            stats = regionprops(this.L, 'Area');
        else
            disp('implement version without label matrix');
        end
        
        for ci = 1:length(this.bonds{t})

            l = this.bonds{t}(ci).label;
            
            geom = struct();

            % the region props ones
            if ~isempty(this.L) 
                
                if ~isnan(l)
                    % bonds are single pixel thickness, so area = length
                    geom.length = stats(l).Area;
                else
                    geom.length = NaN;
                end
            end
            
            % set property
            this.bonds{t}(ci).setGeometry(geom);
        end
            
    %----------------------------
	% curved case
    %----------------------------
    
    else
        
        error('curved case not implemented');
    end
end
