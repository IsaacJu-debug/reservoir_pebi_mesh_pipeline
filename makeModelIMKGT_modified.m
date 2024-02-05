%% Faulted 2.5D reservoir grid
% In this example, we use the upr module to construct a 2.5 D faulted
% reservoir grid.
clear all;
close all;
clc;

plottingFlag = true;

mrstModule add upr

%% Plotting functionality
fig2D = @() figure('Position', [0,0,800,500]);
fig3D = @() figure('Position', [0,0,1300, 650]);
alpha = 0.6;
cmap  = jet*alpha + (1-alpha);
setAxProps = @(ax) set(ax, 'View'              , [-10, 45]         , ...
                           'PlotBoxAspectRatio', [4.40, 1.86, 1.00], ...
                           'Projection'        , 'Perspective'     , ...
                           'Box'               , 'on'              , ...
                           'ColorMap'          , cmap              );

%% Load constraining curves
% We start by loading a data structure that contains points describing the
% outline, faults, and well positions
pth = fullfile(mrstPath('upr'), 'datasets');
rp  = load(fullfile(pth, 'reservoirPoints.mat'));
rp  = rp.reservoirPoints;

xMin = min( rp.outline(:,1) );
xMax = max( rp.outline(:,1) );
yMin = min( rp.outline(:,2) );
yMax = max( rp.outline(:,2) );

rp.outline = [ xMin, yMin; ...
               150, yMin; ...
               533, yMin; ...
               733, yMin; ...
               xMax, yMin; ...
               xMax, yMax; ...
               870, yMax; ... 
               670, yMax; ... 
               425, yMax; ... 
               xMin, yMax ];

rp.faultLines{1} = [ 870, yMax; 
                     rp.faultLines{1}; ...
                     150, yMin];
rp.faultLines{2} = [ 733, yMin; 
                     rp.faultLines{2}; ...
                     670, yMax];
rp.faultLines{3} = [ 425, yMax; 
                     rp.faultLines{3}; ...
                     533, yMin];      
     
if plottingFlag                   
  fig2D(), hold on
  plot(rp.outline(:,1), rp.outline(:,2), 'ko')
  plotLinePath(rp.faultLines,'b-+');
  plotLinePath(rp.wellLines,'.r', 'markerSize', 20);
  box on, axis equal tight
end 

%% Generate 2D PEBI grid
% We construct a 2D PEBI grid from the points using pebiGrid, with
% refinement around the wells
rng(2019)
n   = 25; % Approximate number of cells in x-direction
L   = max(rp.outline);
G2D = pebiGrid2D(max(L)/n, L          , ...
    'polybdr'        , rp.outline   , ... % Outline
    'faceConstraints', rp.faultLines, ... % Fault lines
    'FCFactor'       , 0.8          , ... % Relative size of fault cells
    'cellConstraints', rp.wellLines , ... % Well coordinates
    'CCRefinement'   , true         , ... % Refine
    'CCFactor'       , 0.1          , ... % Relative size of well cells
    'interpolateFC'  , true         , ... % Interpolate along fault lines
    'CCEps'          , 0.08*max(L)  );    % Refinement transition
G2D = removeShortEdges(G2D, 1); % The grid may contain very short edges.
                                % We remove these
                                
%% Plot 2D grid
% Well cells are identified with G.cells.tag = true. We will need the well
% number later, so we find and store these.
wellNo2D                = nan(G2D.cells.num,1);
wellNo2D(G2D.cells.tag) = 1:nnz(G2D.cells.tag);

if plottingFlag
  fig2D(), plotGrid(G2D); axis equal tight, box on                 % Grid
  plotFaces(G2D, G2D.faces.tag, 'edgeColor', 'b', 'lineWidth', 2); % Faults
  plotGrid(G2D, G2D.cells.tag , 'faceColor', 'r');                 % Wells
end 
  
%% Find regions
% The faults divide the reservoir into six distinct regions. We identify
% these using functionality from the coarsegrid module
mrstModule add coarsegrid
p = ones(G2D.cells.num,1);
p = processPartition(G2D, p, find(G2D.faces.tag));

%% Make 2.5D reservoir model
% We construct a volumetric reservoir model by extruding the 2D grid using
% makeLayeredGrid
nLayers = 23; % Number of layers
layerThickness = ones(nLayers,1)*5 + (rand(nLayers,1)*2 - 1)*2;
G0           = makeLayeredGrid(G2D, layerThickness);
G0           = computeGeometry(G0);
G0.cells.tag = repmat(G2D.cells.tag, nLayers, 1);
G0.faces.tag = false(G0.faces.num,1);
G0.faces.tag(abs(G0.faces.normals(:,3))<.01) = repmat(G2D.faces.tag, nLayers, 1);
wellNo       = repmat(wellNo2D, nLayers, 1);
layerID      = reshape(repmat(1:nLayers, G2D.cells.num, 1), [], 1);
compartID    = repmat(p,nLayers,1);

%% Plot the resulting layered grid
prm = randperm(nLayers)';

if plottingFlag
  fig3D(); plotCellData(G0, prm(layerID), 'edgealpha', 0.2)
  outlineCoarseGrid(G0, compartID,'EdgeColor','w','LineWidth',2);
  setAxProps(gca), camlight();
  colormap(jet)
  axis off
end 
  
%% Remove cells
% To add more realism to the reservoir, we mimic erosion and inaccessible
% parts of the formation by removing cells
bnds = [0 8; 1 7; 1 7; 5 3; 2 6; 15 3];
flag = false(G0.cells.num,1);
for i=1:6
    flag = flag | ((compartID==i) & ...
        (layerID<=bnds(i,1) | layerID>=(nLayers-bnds(i,2))));
end
[G, cellMap] = removeCells(G0, flag);
G      = computeGeometry(G);
layer  = layerID(cellMap);
parts  = compartID(cellMap);
wellNo = wellNo(cellMap);

%% Populate with petrophysical properties
% We populate the grid with petrophysical properties drawn from a layered,
% lognormal, isotropic permeability field. To fake faults with
% displacement, the permeability inside each compartment is sampled from
% the same cube. This produces a plausible layering structure, but does not
% preserve the areal correlation within a single geological layer on
% opposite sides of a fault. 
permMean = [10, 912, 790, 90, 10];   % Mean permeability in each layer
N        = [90, 30, nLayers];        % Num of points in each axial direcion
ind      = [1,5,13,15,20,nLayers+1]; % Layer indices
K        = reshape(logNormLayers(N, permMean, 'indices', ind), N);

perm  = nan(G.cells.num,1);
for i = 1:6
    idx = parts==i;
    perm(idx) = sampleFromBox(G, K, find(idx))*milli*darcy;
end
lperm = log10(perm);
poro = (lperm - min(lperm))./(max(lperm) - min(lperm)).*0.7 + 0.1;
rock = makeRock(G, perm, poro);

%% Plot permeability
if plottingFlag
  fig3D(); plotCellData(G, log10(rock.perm), 'edgealpha', 0.2)
  setAxProps(gca), camlight();
end

%% Shift vertical coordinates
% Finally, we shift the vertical coordinates of the grid to mimc geological
% activity.
if true
  x    = G.nodes.coords;
  xmax = max(x);
  xmin = min(x);
  xr   = (x - mean(x))./((xmax - xmin)/2)*2; % Transform to reference domain
  z    = peaks(xr(:,1), xr(:,2));            % Vertical shift
  x(:,3) = x(:,3) - mean(x(:,3));
  % We shift the vertical coordinates by using the built-in MATLAB function
  % peaks, and also make it thicker along the y axis
  x(:,3) = (x(:,3) + z*7).*(x(:,2)./xmax(2)+1).^2*0.25;
  x(:,3) = x(:,3) - min(x(:,3)) + 1000*meter; % Normal reservoir depth
  G.nodes.coords = x;
  G = computeGeometry(G);
end
  
%% Make wells
W = [];
for wNo = 1:max(wellNo)
    cells = find(wellNo == wNo);
    W     = addWell(W, G, rock, cells, 'name', '');
end
if plottingFlag
  pw = @() plotWell(G, W, 'color', 'k', 'height', 25, 'LineWidth', 3);
end 

%% Plot the final result
% View 1: with wells
if plottingFlag
  fig3D(), plotCellData(G, log10(rock.perm), 'edgealpha', 0.2); pw();
  setAxProps(gca), camlight, view([-10,45]);
end 

% View 2: exploded view of compartments
if plottingFlag
fig3D()
xmid = mean(G.cells.centroids);
  for i = 1:6
      g = G;
      xlmid = mean(G.cells.centroids(parts==i,:));
      v = xlmid - xmid;
      g.nodes.coords = g.nodes.coords + 0.3*v;
      plotCellData(g, log10(rock.perm), parts==i, 'edgealpha', 0.2);
  end
  setAxProps(gca); camlight
  axis tight
end

% Allocating space for polyhedra output
numFacesPerCell = diff(G.cells.facePos);
numNodesPerFace = diff(G.faces.nodePos);

nPrisms = zeros(9,1 ); %prism3 (aka wedge), ..., prism11
for nSides = 3:11
  nPrisms(nSides-2) = sum(numFacesPerCell == ( nSides + 2) );
end

if sum(nPrisms) ~= G.cells.num
  error("Cells are not just prism3 to prism11");
end

storage = 0;
for nSides = 3:11
  nFaces = nSides + 2;
  offset = nFaces + 1 + 2*nSides + (nFaces-2)*4;
  storage = storage + nPrisms(nSides-2)*offset;
end

% Populate connectivity array
connectivity = zeros(storage+G.cells.num,1);

ptr = 1;
for iCell = 1:G.cells.num
  
  nFaces = numFacesPerCell(iCell);
  nSides = nFaces - 2;
  offset = nFaces + 1 + 2*nSides + (nFaces-2)*4;
  
  connectivity( ptr ) = numFacesPerCell(iCell);
  
  tmp = [ offset; nFaces ];
  
  cellFaces = G.cells.faces(G.cells.facePos(iCell) : G.cells.facePos(iCell+1)-1, 1);
  
  for iFace = 1:nFaces
    faceID = cellFaces(iFace);
    N = numNodesPerFace( faceID );
    nodes = G.faces.nodes(G.faces.nodePos(faceID) : G.faces.nodePos(faceID+1)-1, :) - 1;

    
    tmp = [tmp; N; nodes];

  end
  
  connectivity(ptr:ptr+offset) = tmp;
  
  ptr = ptr + offset+1;  

end
  
%%%%%% 
% OUTPUT vtk
% Open output file
FID_vtk = fopen('pebi_mrst.vtk','w');

% Header
fprintf(FID_vtk,'%s\n','# vtk DataFile Version 4.2');
fprintf(FID_vtk,'%s\n','pebi grid mesh');
fprintf(FID_vtk,'%s\n','ASCII');
fprintf(FID_vtk,'%s\n','DATASET UNSTRUCTURED_GRID');

fprintf(FID_vtk,'POINTS %d float\n',G.nodes.num);
fprintf(FID_vtk, ...
        '%13.7e %13.7e %13.7e\n', ...
        G.nodes.coords' );

fprintf(FID_vtk,'CELLS %d %d\n',G.cells.num, length(connectivity) );
fprintf(FID_vtk, ...
        '%d %d %d %d %d %d %d %d %d %d\n', ...
        connectivity );      
fprintf(FID_vtk, '\n');      
      
      
fprintf(FID_vtk,'CELL_TYPES %d\n',G.cells.num);
fprintf(FID_vtk, ...
        '%d\n', ...
        42*ones(G.cells.num,1) );      
 fprintf(FID_vtk, '\n');      
      
      
fprintf(FID_vtk,'CELL_DATA %d\n',G.cells.num);
fprintf(FID_vtk,'SCALARS PERM float 1\n');
fprintf(FID_vtk,'LOOKUP_TABLE default\n');
fprintf(FID_vtk, ...
        '%11.5e\n', ...
        rock.perm );      
fprintf(FID_vtk,'SCALARS PORO float 1\n');
fprintf(FID_vtk,'LOOKUP_TABLE default\n');
fprintf(FID_vtk, ...
        '%11.5e\n', ...
        rock.poro ); 
      
fclose( FID_vtk);

