%% 2.5D reservoir grid
clear all;
close all;
clc;

mrstModule add upr


%% Load domain outline

outline   = load('outline.xy');
wellLines = load('wellLines.xy');
nCases = size(wellLines,1);

nLayers = 256;%1; % Number of layers
layerThickness = ones(nLayers,1)*25 + (rand(nLayers,1)*2 - 1)*2;

%% permeability field
rng(2019,'twister')
prm = randperm(nLayers)';
permMean = [1 ];                 % Mean permeability in each layer
N        = [100 , 100, nLayers]; % Num of points in each axial direcion
ind      = [1,nLayers+1];        % Layer indices
K        = reshape(logNormLayers(N, permMean, 'indices', ind), N);

for iCase = 1:nCases
  
  disp(strcat('Processing case # ', num2str(iCase,'%0.4d'),' / ',num2str(nCases,'%0.4d')));
  
  %% Generate 2D PEBI grid
  % 2D PEBI grid is constructed for the provided domain outline, with
  % refinement around the wells
  n   = 110; % Approximate number of cells in x-direction
  domainSize   = [ max(outline(:,1)) - min(outline(:,1)), ...
    max(outline(:,2)) - min(outline(:,2)) ];
  
  [G2D, Pts ] = pebiGrid2D( ...
    max(domainSize)/n, ...
    domainSize          , ...
    'polybdr'        , outline , ...
    'cellConstraints', wellLines(iCase,:) , ...
    'CCRefinement'   , true         , ...
    'CCFactor'       , 0.1          , ...
    'CCEps'          , 0.08*max(domainSize)  );
  G2D = removeShortEdges(G2D, 1);
  
  % Well cells are identified with G.cells.tag = true.
  wellNo2D                = nan(G2D.cells.num,1);
  wellNo2D(G2D.cells.tag) = 1:nnz(G2D.cells.tag);
  
  
  %% Make 2.5D reservoir model
  G           = makeLayeredGrid(G2D, layerThickness);
  G           = computeGeometry(G);
  G.cells.tag = repmat(G2D.cells.tag, nLayers, 1);
  G.faces.tag = false(G.faces.num,1);
  wellNo      = repmat(wellNo2D, nLayers, 1);
  
  %% Populate with petrophysical properties
  % We populate the grid with petrophysical properties drawn from a layered,
  % lognormal, isotropic permeability field.
  perm = 100*sampleFromBox(G, K, 1:G.cells.num )*milli*darcy;
  lperm = log10(perm);
  poro = (lperm - min(lperm))./(max(lperm) - min(lperm)).*0.3 + 0.1;
  rock = makeRock(G, perm, poro);
  
  min(perm)
  max(perm)
  min(poro)
  max(poro)

  size(poro)

  %% Write output
  %wellCellID = write_well_box( iCase, G, wellNo );    
  %write_well_index( iCase, G, rock, wellCellID );  
  write_vtk( G, Pts, rock, wellNo, strcat( num2str(iCase,'%0.4d'), '_domain.vtk' ));
  
end
