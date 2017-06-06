function [VoronoiLattice,vertices] = ExtractLattice(memseg, varargin)
    % EXTRACTLATTICE Extract a Voronoi structure from a membrane segmentation
    % 
    % [VoronoiLattice, vertices] = ExtractLattice(memseg)
    % [VoronoiLattice, vertices] = ExtractLattice(memseg, options)
    %
    % memseg:           image containing membrane segmentation 
    %                   should have membrane 1, cells 0, and, if applicable
    %                   background -1 
    % options:          struct containing fields:
    % - closeSize:      size of disk with which to close holes
    % - areaOpenSize    minimal size of cells
    % - minVertexDist   merge vertices closer together than this
    % - bgCells         return disconnected background areas as cells
    % 
    % VoronoiLattice:   struct containing
    % - vertexPosition: vertices [x y] 
    % - cellVertices:   counterclockwise ordered bonds 
    % - bulkIdx:        indices of bulk cells
    % - bdryIdx:        indices of bdry cells
    % - bgIdx:          indices of background "cells"
    % - L:              label matrix of cleaned up segmentations    
    %
    % vertices: image of vertices for debugging
    %
    %---------------------------------------------------------------------
    % Copyright 2014 Idse Heemskerk and Sebastian Streichan
    %---------------------------------------------------------------------
    
    %-----------------------------
    % initialize options structure
    %-----------------------------
    
    if nargin==2
        options = varargin{1};
    else
        options = struct();
    end
    
    if isfield(options, 'closeSize')
        closeSize = options.closeSize;
    else
        closeSize = 5;
    end
    if isfield(options, 'areaOpenSize')
        areaOpenSize = options.areaOpenSize;
    else
        areaOpenSize = 10;
    end
    if isfield(options, 'minVertexDist') 
        r = options.minVertexDist;
    else
        r = 0;
    end    
    if isfield(options, 'bgCells') 
        bgCells = options.bgCells;
    else
        bgCells = true;
    end    
    
    %-----------------------------
    % clean up segmentation
    %-----------------------------
    
    membrane = imclose(memseg, strel('disk',closeSize));
    
    % for a binary image this would do:
    % membrane = ~bwareaopen(~membrane, areaOpenSize, 4);   
    % but for an image with bg = -1, we need to not lose the distinction
    % between background, membrane and cell
    % ~membrane is equivalent to membrane~=1 so we get bg & cell
    removed = ~membrane & ~bwareaopen(~membrane, areaOpenSize, 4);
    membrane(removed) = 1;

    % from here on membrane is binary and we make a separate image of
    % background 
    background = membrane == -1;
    membrane = membrane == 1;
    
    % 'thin' doesn't remove edges that
    % cross the boundary, except when and edge is "tangent"
    % thick frame around the image so in all case bonds stay connected to boundary
    d=6;
    tmp = true(size(membrane) + d);
    tmp(d/2+1:end-d/2,d/2+1:end-d/2) = membrane;    
    tmp = bwmorph(tmp, 'thin', 'inf');
    
    % shrink removes edges with one vertex only (lines sticking out)
    tmp = bwmorph(tmp, 'shrink', 'inf');
    
    % remove frame 
    membrane = tmp(d/2+1:end-d/2,d/2+1:end-d/2);
    membrane = bwmorph(membrane, 'clean');
    
    % make a binary image of the vertices
    % introduce artificial vertices where edges intersect the boundary
    % boundary vertices may end up with a one pixel displacement but their
    % location is not crucial anyway, they're a bookkeeping device only
    bdrymembrane = membrane;
    bdrymembrane(:,1)=1;
    bdrymembrane(1,:)=1;
    bdrymembrane(end,:)=1;
    bdrymembrane(:,end)=1;
    bdrymembrane = bwmorph(bdrymembrane, 'thin', 'inf');
    
    vertices = bwmorph(bdrymembrane, 'branchpoints');
    
    % make binary images of cells, separate into cells intersecting the
    % image boundary and interior cells
    allcells = ~membrane;
    % make a label matrix for all cells
    cellsL = bwlabel(allcells,4);
    % identify bg "cells"
    bgIdx = unique(background.*cellsL);
    
    % if bg cells are not wanted, remove them
    if ~bgCells
        for i = setdiff(bgIdx',0)
            cellsL(cellsL == i) = 0;
        end
        allcells = cellsL>0;
        cellsL = bwlabel(allcells,4);
        %bulkIdx = setdiff(bulkIdx, bgIdx);
        bgIdx = 0;
    end

    nCells = max(cellsL(:));
    
    % separate labels into bdry and bulk
    bulkcells = imclearborder(allcells, 4);
    bdrycells = allcells & ~bulkcells;
    bdryIdx = setdiff(unique(cellsL.*bdrycells),0);
    bulkIdx = setdiff(unique(cellsL.*bulkcells),0);
    
    %------------------------------
    % create the Voronoi structure
    %------------------------------
    
    % get the vertex coordinates
    [vertexY, vertexX] = ind2sub(size(vertices), find(vertices));
    V = [vertexX, vertexY];

    % make a cell array tracking which vertices belong to which cells
    C = cell([1 nCells]);

    % for each vertex, add vertex index to cells that are 8-connected to it
    ysize = size(cellsL,1);
    xsize = size(cellsL,2);
    
    for i=1:length(vertexX)
        
        belongsto = unique(cellsL(max(vertexY(i)-1,1):min(vertexY(i)+1,ysize),...
                             max(vertexX(i)-1,1):min(vertexX(i)+1,xsize)));
        for vertj = 2:length(belongsto)
            C{belongsto(vertj)} = [C{belongsto(vertj)}, i];
        end
    end

    % Make sure cells are sorted!
    % PROBLEM: this only works for convex cells
    for I = 1 : nCells
        nVertices = length(C{I});
        [~, Index] = ClockWiseSortCell( C{I}(1:nVertices), V);
        %cellVertices{I} = [cellVertices{I}(Index), cellVertices{I}(Index(1))];
        C{I} = C{I}(Index);
    end
    
    %---------------------------------------------------
    % merge vertices that are close together
    %---------------------------------------------------
    
    % with r=0 this removes vertices that are adjacent pixels
    
    % make a label matrix of dilated vertices
    dilVerts = imdilate(vertices, strel('disk',r));
    vertL = bwlabel(dilVerts);

    % make an image of the interfaces
    % add the interfaces to the cell labelmatrix
    if r==0
        bonds = membrane & ~imdilate(vertices,strel('square',3));
    else
        bonds = membrane & ~dilVerts;
    end
    bondsL = bwlabel(bonds);
    bondsL = uint16(bondsL) + nCells.*uint16(bondsL>0);
    %bondsL = imdilate(bondsL,strel('square',3));
    bondIdx = (1:max(bondsL(:))) + nCells;
	totalL = uint16(cellsL);
    totalL(bondsL>0) = bondsL(bondsL>0);
    
    % a vector of mappings from vertices onto dilated (and potentially
    % merged) vertices
    mergeIdx = vertL(sub2ind(size(vertL), V(:,2), V(:,1)));

    % positions of merged vertices
    stats = regionprops(vertL, 'Centroid');
    redV = round(cat(1,stats(:).Centroid));

    % get bond labels attached to vertex
    % LOT OF POTENTIAL TO MAKE THIS FASTER
    dilLVerts = imdilate(vertL,strel('square',3));
    vertBondLabels = {};
    for i = 1:max(vertL(:))
        vbondL = unique(bondsL(dilLVerts == i));
        vbondL(vbondL==0)=[];
        vertBondLabels{i} = vbondL;
    end
    
    % cells referencing the reduced vertex list
    redC = {};
    for i = 1:length(C)
        redC{i} = mergeIdx(C{i});
    end

    C = redC;
    V = redV;
    
    % add third column for generality
    V = [V, 0*V(:,1)];
    
	% after merging the same vertex can appear twice in the cell, so we
	% remove doubles
    for i = 1:nCells
        C{i} = unique(C{i},'stable');
    end
    
    %---------------------------------------------------

    VoronoiLattice = struct('vertexPosition', {V},...
                            'cellVertices', {C},...
                            'bdryIdx', bdryIdx,...
                            'bulkIdx', bulkIdx,...
                            'bgIdx', bgIdx,...
                            'bondIdx', bondIdx,...
                            'nCells', nCells,...
                            'vertBondLabels', {vertBondLabels},...
                            'L', totalL);
end