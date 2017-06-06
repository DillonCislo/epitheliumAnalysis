%------------------------------------------------------
% script to test lattmin conversion
%------------------------------------------------------

close all;
clear all;

[testScriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(testScriptPath);
cd('..');
addpath(genpath(pwd));

%----------------------------------------------
% create Voronoi structure from data 
%----------------------------------------------

image = imread('zyx17A.actGFP_series.2A_Ecad_MIP.tif');
memseg = ~imread('zyx17A.actGFP_series.2A_Ecad_seg.tif');

options = struct('closeSize', 0, 'minVertexDist', 0, 'areaOpenSize', 50);
[VoronoiLattice,vertices] = ExtractLattice(memseg, options);

% inspect Voronoi lattice

imshow(image,[]);

V = VoronoiLattice.vertexPosition;
C = VoronoiLattice.cellVertices;

cellIdx = VoronoiLattice.bulkIdx';

for i = cellIdx
    patch(V(C{i},1), V(C{i},2), [1 0 0], 'FaceColor','none','EdgeColor','green','LineWidth',1)
end

%%
%----------------------------------------------
% create Latmin structure from Voronoi struct
%----------------------------------------------

bulkCellIdx = setdiff(VoronoiLattice.bulkIdx, VoronoiLattice.bgIdx);
trim = true;
g = GLattConversion2(C(bulkCellIdx), V, trim);

options = struct('cellIndex', true);
figure, imshow(image);
LatticePresentation2(g, options)

%%
%----------------------------------------------
% create cellLayer from lattmin
%----------------------------------------------

% create a CellLayer object called cellLayer
cellLayer = CellLayer();

t = 1;
cellLayer.initTime(t, 'Lattmin', g);
cellLayer.setCellLabels(t, bulkCellIdx(g.cti2ci));
cellLayer.setLabelMatrix(VoronoiLattice.L);

%%
% check that conversion back and forth didn't destroy anything 
g2 = cellLayer.getLattmin(t);
LatticePresentation2(g2, options)

%%
% check that some methods work
mycell = cellLayer.getCell(t, 15)
cellNN = mycell.getNeighbors();
cellNNlabels = [cellNN.label]
