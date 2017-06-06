function initTime(this, t, inputType, input, options)
    % INITTIME initialize the cell layer at some particular time
    %
    % initTime(inputType, input)
    %
    % inputType:    string 'Lattmin', 'Kevin' or 'image'
    % input:        data on the format specified by inputType

    assert(isnumeric(t), 'first argument should be time');
    
    if strcmp(inputType, 'Lattmin')

        initTimeFromLattmin(this, t, input);

    elseif strcmp (inputType, 'Kevin')

        initTimeFromKevin(this, t, input(t));

    elseif strcmp (inputType, 'image')

        initTimeFromImage(this, t, input, options);
        
    else
        error(['input type ' inputType ' not recognized, try: Lattmin, Kevin, image']);
    end
    
    % initialize the state to NaN
    if ~isempty(this.cellState)
        for i = 1:numel(this.cells{t})
            this.cells{t}(i).setState(nan([1 numel(this.cellState)]));
        end
    end
    if ~isempty(this.bondState)
        for i = 1:numel(this.bonds{t})
            this.bonds{t}(i).setState(nan([1 numel(this.bondState)]));
        end
    end
end

function initTimeFromLattmin(this, t, g)
    % INITTIMEFROMLATTMIN init at some time from lattmin structure
    %
    % initTimeFromLattmin(t, g)
    %
    % t:    time
    % g:    Lattmin structure, which has field:
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

    disp('initializing from lattmin');

    % cells
    %---------------

    nCells = length(g.cells);
    CeLLArray(nCells) = CeLL; % CeLL object array at time t

    for i = 1:nCells
        % CeLL(cellLayer, time, label, bondInd)
        CeLLArray(i) = CeLL(this, t, i, g.cells{i});
    end

    this.cells{t} = CeLLArray;

    % bonds
    %---------------

    nBonds = size(g.bonds,1);
    BondArray(nBonds) = Bond; % Bond object array at time t

    bondVerts = g.bonds(:,1:2);
    bondCells = g.bonds(:,3:4);

    for i = 1:nBonds
        % Bond(cellLayer, time, label, vertInd, cellInd)
        BondArray(i) = Bond(this, t, i, bondVerts(i,:), bondCells(i,:));
    end

    this.bonds{t} = BondArray;

    % vertices
    %---------------

    if size(g.verts,2) == 2
        this.vertices2{t} = g.verts;
    else
        this.vertices{t} = g.verts;
    end
end

%-----------------------------------------------------------------

function initTimeFromImage(this, t, seg, options)
    % INITTIMEFROMIMAGE init at some time from segmented image
    % 
    % initTimeFromImage(t, seg, options)
    %
    % seg:          image containing membrane segmentation 
    %               should have membrane 1, cells 0, and, if applicable
    %               background -1 
    %
    % options:      struct containing fields:
    %
    % closeSize:    size of disk with which to close holes
    % areaOpenSize  minimal size of cells
    % minVertexDist merge vertices closer together than this
    % bgCells       return disconnected background areas as cells
    % trim:         trim lattice to create boundary bonds, this removes
    %               outside cells with twofold vertices and only keeps the
    %               edges attached to other cells, useful for mech inverse

    VoronoiLattice = ExtractLattice(seg, options);

    V = VoronoiLattice.vertexPosition;
    C = VoronoiLattice.cellVertices;

    if isfield(options, 'trim')
        g = GLattConversion2(C, V, options.trim);
    else
        g = GLattConversion2(C, V, false);
    end

    initTimeFromLattmin(this, t, g);
    
    this.setCellLabels(t, g.cti2ci); 
    this.setLabelMatrix(t, VoronoiLattice.L)

    % determine bond labels
    % bonds are oriented, labels are not, so two opposing bonds have the
    % same label

    disp('determining bond labels');

    nBonds = length(this.bonds{t});
    cnt = 0;

    for i = 1:nBonds

        bVerts = this.bonds{t}(i).vertInd;
        bLabel = intersect(VoronoiLattice.vertBondLabels{bVerts});

        if isempty(bLabel)
            % happens on fake bonds on the image edge

            bLabel = NaN;

        elseif numel(bLabel) == 2
            % happens for two sided cells

            cnt = cnt + 1;
            dilBdry = imdilate(this.L == bLabel(1),strel('square',3));
            touchingCellLabels = setdiff(this.L(dilBdry), bLabel(1));

            if numel(touchingCellLabels) == 2
                if any(this.bonds{t}(i).cellInd == 0)
                    bLabel = bLabel(1);
                else
                    bLabel = bLabel(2);
                end
            else
                if any(this.bonds{t}(i).cellInd == 0)
                    bLabel = bLabel(2);
                else
                    bLabel = bLabel(1);
                end
            end

        elseif numel(bLabel) > 2
            % happens for pair of two sided cells

            disp('work out the bond label of pair of two sided cells!');
            bLabel = NaN;
        end
        this.bonds{t}(i).setLabel(bLabel);
    end
end

%-----------------------------------------------------------------

function initTimeFromKevin(this, t, g)
    % INITTIMEFROMKEVIN init at some time from kevin's structure
    %
    % initTimeFromKevin(t, g)
    %
    % t:    time
    % g:    kevin's structure

    disp('initializing from kevin');

    nCells = length(g.Cdat);
    nVerts = length(g.Vdat);
    nBonds = nCells + nVerts; 

    % Build needed structures to initialize verts, cells, & bonds.
    bondVerts = zeros(nBonds,2);
    bondCells = zeros(nBonds,2);
    this.vertices2{t} = zeros(nVerts,3);
    nb = 1;
    for v = 1:nVerts
        this.vertices2{t}(v,1:2) = [g.Vdat(v).vertxcoord, g.Vdat(v).vertycoord];
        for nv = g.Vdat(v).nverts
            if ( ~sum( (bondVerts(:,1) == nv) .* (bondVerts(:,2) == v) ) )
                bondVerts(nb,1) = v;
                bondVerts(nb,2) = nv;
                bondCells(nb,:) = g.Vdat(v).ncells(ismembc(g.Vdat(v).ncells,g.Vdat(nv).ncells));
                nb  = nb + 1;
            end
        end
    end

    if (nb <= nBonds)
        bondVerts(nb:nBonds,:) = [];
        bondCells(nb:nBonds,:) = [];
        nBonds = size(bondCells,1);
    end

    % CeLL object array at time t
    CeLLArray(nCells) = CeLL;

    for i = 1:nCells
        cellBonds = zeros(1,length(g.Cdat(i).orderednverts)-1);
        for j = 1:length(cellBonds)
            v = sort([g.Cdat(i).orderednverts(j),g.Cdat(i).orderednverts(j+1)]);
            cellBonds(j) = find( (bondVerts(:,1) == v(1)) .* (bondVerts(:,2) == v(2)) );
        end
        CeLLArray(i) = CeLL(this, t, i, cellBonds);
    end

    this.cells{t} = CeLLArray;

    % Bond object array at time t
    BondArray(nBonds) = Bond;

    for i = 1:nBonds
        BondArray(i) = Bond(this, t, i, bondVerts(i,:), bondCells(i,:));
    end

    this.bonds{t} = BondArray;
end