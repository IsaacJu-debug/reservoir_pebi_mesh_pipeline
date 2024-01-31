%% 2.5D reservoir grid
%% This script is to generate fractured reservoir models with heterogeneous geophysical properties
%% Currently, the geomodel is given as an input created from a python file.

clear all;
close all;
clc;

mrstModule add upr
%% flags for reading files
isFractured = true;
isSplited = true;
isDebug = false;

filePad = 'simple';
filename = 'pebiMultiphasePoromechanics.xml';
baseFileName = 'pebiMultiphasePoromechanics_base.xml';
%% Load domain outline
outline   = load('outline.xy');
wellLines = load('wellLines.xy');

if isFractured
  faceLines = {[26000.0,23000.0;29100.0,29200.0], [11600.0,25200.0;25600.0,43300.0]}; % predefined fault lines
  
  if isSplited
    fprintf('Generate fractures along fault lines\n');
    % need to use thomas' mesh_doctor to generate fracture
    filePad = 'splited';
  else
    % this is for permeability multiplier
    fprintf('Use permeability multiplier for mimicing seal faults\n');
    filePad = 'frac';
  end
end
kMapSizes = 30; % the grid resolution of geomodel

reservoirTop = 2700;
reservoirBottom = 2500;
overburdenPerm = 1e-20;
underburdenPerm = 1e-20; 
overburdenPoro = 0.15; % avoid too small values
underburdenPoro = 0.15; % avoid too small valuees 

coef = 1.0; % 1, 2.1, 4.4, 8.9
amplificationFactor = 1.85;

nCases  = 1;
Lz = 5200;

nLayersInReservoir = 10*coef;
nLayersInUnderburden = 3*coef;
nLayersInOverburden = 3*coef;

dzReservoir = (reservoirTop-reservoirBottom) / nLayersInReservoir;

done = false;
dz = dzReservoir;
top = reservoirTop + dz;
dzOverburden = [dz];
maxDzOverburden = (Lz-reservoirTop) / nLayersInOverburden;
while ~done
  dz = amplificationFactor * dz
  if dz > maxDzOverburden
    dz = maxDzOverburden;
  end
  if top + dz > Lz
    dz = Lz - top;
    if dzOverburden(length(dzOverburden)) > dz
      dzOverburden(length(dzOverburden)) = dzOverburden(length(dzOverburden)) + dz; 
    else
      dzOverburden = [dzOverburden,dz];
    end
    break;
  end
  dzOverburden = [dzOverburden,dz];
  top = top + dz;
end
dzUnderburden = dzOverburden([length(dzOverburden):-1:1]);    
nLayersInOverburden = size(dzOverburden,2);
nLayersInUnderburden = size(dzUnderburden,2);

dzUnderburden
dzOverburden

nLayers = nLayersInReservoir + nLayersInUnderburden + nLayersInOverburden; % Number of layers
%nFracLayers = ;
%nWellLayers = ;
layerThickness = ones(nLayers,1);
layerThickness([1:nLayersInUnderburden]) = dzUnderburden(1,1:nLayersInUnderburden);
layerThickness([(nLayersInUnderburden+1):(nLayersInUnderburden+nLayersInReservoir)]) = ones(nLayersInReservoir,1)*dzReservoir;
layerThickness([(nLayersInUnderburden+nLayersInReservoir+1):nLayers]) = dzOverburden(1,1:nLayersInOverburden);

if isDebug
    fprintf('Debugging mode. Change the layer number to 1\n');
    layerThickness = ones(1,1);
    layerThickness(1) = 2;
    filePad = 'debug_splited';
end

layerThickness

%% permeability field
rng(2019,'twister')
logk_mean = log(0.1); % k means 0.1 mD 
logk_std = 2.5;
nX = kMapSizes;
nY = nX;
nZ = 260;

N   = [25*coef, 25*coef, nLayers]; % Num of points in each axial direcion
ind = [1,nLayers+1]; % Layer indices

%ax = randi([floor(1*nX/2),floor(3*nX/4)], 1, 1); % correlation length is a discrete uniform distribution
%ay = randi([floor(1*nY/2),floor(3*nY/4)], 1, 1); % the isotropic ratio between ax and az is 1 in this setting           
%kmap = vonk2d_fixed( 1, 1, 1, ax(1), ay(1), nX, nY, 1, 3, 1, (1), (1) );
%K = reshape(kmap(:, :), [nX , nY, 1] );
%K = scaleLogK(K, logk_mean, logk_std);

K = dlmread('permeability30x30x1000.out');
K = reshape(K, [nX , nY, nZ]);

savePcolor(log10(K(:,:,1)), './', 'permeabilityMap1.png' );
savePcolor(log10(K(:,:,2)), './', 'permeabilityMap2.png' );

for iCase = 1:nCases
  
  %% Generate 2D PEBI grid
  % 2D PEBI grid is constructed for the provided domain outline, with
  % refinement around the wells
  n = N(1); % Approximate number of cells in x-direction
  domainSize = [ max(outline(:,1)) - min(outline(:,1)), ...
    max(outline(:,2)) - min(outline(:,2)) ];
  
  if isFractured
    % Special operations when isFractured is true
    [G2D, Pts] = pebiGrid2D( ...
        max(domainSize)/n, ...
        domainSize          , ...
        'polybdr'        , outline , ...
        'cellConstraints', wellLines ,...
        'CCRefinement'   , true         , ...
        'CCFactor'       , 0.2          , ...
        'CCEps'          , 0.08*max(domainSize), ...
        'faceConstraints', faceLines, ...
        'useMrstPebi', true);
  else
    % Original code for non-fractured scenario
    [G2D, Pts] = pebiGrid2D( ...
        max(domainSize)/n, ...
        domainSize       , ...
        'polybdr'        , outline , ...
        'cellConstraints', wellLines , ...
        'CCRefinement'   , true      , ...
        'CCFactor'       , 0.2       , ...
        'CCEps'          , 0.08*max(domainSize)  );
  end
  G2D = removeShortEdges(G2D, 1);
  
  % Well cells are identified with G.cells.tag = true.
  wellNo2D                = nan(G2D.cells.num,1);
  wellNo2D(G2D.cells.tag) = 1:nnz(G2D.cells.tag);

  %% Make 2.5D reservoir model
  G           = makeLayeredGrid(G2D, layerThickness);
  G           = computeGeometry(G);
  %G.cells.tag = repmat(G2D.cells.tag, nLayers, 1);
  %G.faces.tag = repmat(G2D.faces.tag, nLayers, 1); % this can control which layer fractures will be located
  G.faces.tag = false(G.faces.num,1); 
  wellNo      = repmat(wellNo2D, 1, 1); % only handle one cell right now

  %% Populate with petrophysical properties
  % We populate the grid with petrophysical properties drawn from a layered,
  % lognormal, isotropic permeability field.
  perm = sampleFromBox(G, K, 1:G.cells.num );
  lperm = log10(perm);
  poro = (lperm - min(lperm))./(max(lperm) - min(lperm)).*0.2 + 0.1;
  rock = makeRock(G, perm, poro);
  
  min(perm)
  max(perm)
  min(poro)
  max(poro)

  size(poro)

  rock.perm(G.cells.centroids(:,3) > reservoirTop) = overburdenPerm;
  rock.perm(G.cells.centroids(:,3) < reservoirBottom) = underburdenPerm; 
  rock.poro(G.cells.centroids(:,3) > reservoirTop) = overburdenPoro;
  rock.poro(G.cells.centroids(:,3) < reservoirBottom) = underburdenPoro; 

  %% Write output
  wellCellID = write_well_box( iCase, G, wellNo, baseFileName, ...
                                reservoirBottom, reservoirTop);    
  write_well_index( iCase, G, rock, wellCellID ); 

  %% Write output
  currentFolder = pwd;
  mesh_folder = '01_mesh';
  xml_folder = '02_xml';
  rock_folder = '03_rock_prop';

  if ~exist(xml_folder, 'dir')
      mkdir(xml_folder);
      disp(strcat('mkdir xml_folder...'));
  end

  if ~exist(mesh_folder, 'dir')
      mkdir(mesh_folder);
      disp(strcat('mkdir mesh_folder...'));
  end

  if ~exist(rock_folder, 'dir')
      mkdir(rock_folder);
      disp(strcat('mkdir rock_folder...'));
  end
   
  vtk_name = strcat( num2str(iCase,'%0.4d'), '_', filePad, '_domain.vtk' );
  vtk_path = fullfile(currentFolder, mesh_folder, vtk_name);

  if isFractured
    % when generating fractured vtk, we need ensure fractures information
    % to be added both topology and fields.

    %wellIndices = [(nLayersInUnderburden+1):(nLayersInUnderburden+nLayersInReservoir)];
    faultIndices = [(nLayersInUnderburden+1):(nLayersInUnderburden+nLayersInReservoir)]; % go through the reservoir
     % this version do not set well type
    G.cells.tag = repmat(G2D.cells.tag, nLayers, 1); % reservoir layer
    G.faces.tag = repmat(G2D.faces.tag, nLayers, 1);

    [G, G2D] = set_fault_depths(faultIndices, G, G2D, nLayers);
    numFacesAlongFault = length(G.faces.tag(G.faces.tag == 1));
    if ~isSplited
      % save rock properities
      perm = kron(rock.perm',[1. 1. 1. ]');
      cell_rock_prop = cat(1, perm, rock.poro' , G.cells.volumes');

      rock.perm = cat(1, rock.perm, zeros(numFacesAlongFault, 1));
      rock.poro = cat(1, rock.poro, zeros(numFacesAlongFault, 1));

    end
    
    nPrisms = write_vtk_field_data( G, rock, numFacesAlongFault, ...
                                    vtk_path, reservoirBottom, reservoirTop, isSplited);
    updateCellBlocksInXml(filename, nPrisms);

  else
    attribute = rock.poro;
    attribute(:) = 1; 
    attribute( G.cells.centroids(:,3) < reservoirBottom ) = 2;
    attribute( G.cells.centroids(:,3) > reservoirTop ) = 3;
    write_vtk( G, Pts, rock, attribute, wellNo, vtk_path);
  end

  
end
