%------------------------------------------------------
% script to test the template classes
%------------------------------------------------------

close all;
clear all;

[testScriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(testScriptPath));
cd(testScriptPath);

% create a CellLayer object called cellLayer
cellLayer = CellLayer();

% initialize time at the first time with a lattmin structure
g = struct();
g.cells{1} = [1 2 3];
g.bonds = [1 2 1 0; 2 3 1 0; 3 1 1 0];
g.verts = [0.5 0.86; -0.5 0.86; -1 0];

time = 1;

cellLayer.initTime(time, 'Lattmin', g);

cellLayer.getCell(time,1).getNeighbors();