classdef Bond < handle
    % representing a bond in CellLayer

    %---------------------------------------------------------------------
    % properties
    %---------------------------------------------------------------------

    properties (Access = private)

        cellLayer
        time 
    end

    properties (SetAccess = protected)

        label

        vertInd
        cellInd

        state
        geometry        % structure containing geometric properties
                        % such as length and curvature
    end

    %---------------------------------------------------------------------
    % public methods
    %--------------------------------------------------------------------- 

    methods

        %------------------------------------------------------
        % constructor
        %------------------------------------------------------

        function this = Bond(cellLayer, time, label, vertInd, cellInd)
            % initialize Bond
            % 
            % Bond(cellLayer, time, label, vertInd, cellInd)
            %
            % cellLayer:    CellLayer to which the Bond belongs
            % vertInd:      list of vertex indices attached to bond
            % cellInd:      list of cell indices bordering bond

            % to make arrays of an object, it needs to be possible to
            % call the constructor without arguments. matlab sucks
            if nargin==0
                return;
            end

            this.cellLayer = cellLayer;
            this.vertInd = vertInd;
            this.cellInd = cellInd;
            this.time = time;
            this.label = label;
            this.geometry = [];
        end

        %------------------------------------------------------
        % setters
        %------------------------------------------------------

        function setState(this, stateVector)
            % SETSTATE set the state vector
            %
            % setState(stateVector)
            %
            % stateVector length should match stateNames

            if isempty(this.cellLayer.bondState)
                error('state variables were never defined');

            elseif numel(stateVector) ~= numel(this.cellLayer.bondState)

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

            this.geometry = geometry;
        end
    end
end

