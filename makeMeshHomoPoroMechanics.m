%% This script is for generating 2D vertical fractured PEBI mesh (Chap 3)
%% 1. We need to add fault field so that mesh_doctor can identify the plane to split
%% 2. We only generate PEBI mesh, conforming to the fault.
clear all;
close all;
clc;

mrstModule add upr
%% Setup parameters
isNoOutput = true; % if we need to output vpd
useCoreyModel = false; % if we use corey model as relative perm model
readGeomodel = false; % true means read geomodel from saved mat files
poroConst = 0.2; % Gege wen's 2021 value
verboseLevel = 2; % 1. only file name; 2. print out cell number. 
%% Load domain outline
rng(2019,'twister');
kmapSizes = [50]; % the grid resolution of geomodel
nX = kmapSizes(1);
nY = nX;

baseLine = 'vertical_fracture.vtk'; % file name
xmlFilename = 'pebiMultiphasePoromechanics.xml'; % xml input for single inclined fracture
nLayers = 1; % Number of layers
nLayersInUnderburden = 0; % Number of layers belonging to underburden 
nLayersInReservoir = nLayers; % Number of layers belonging to reservoir

thickness = 500.0; % plain deformation
layerThickness = ones(nLayers,1)*thickness/nLayers; % two-body system

reservoirBottom = 0.0;     % the bottomo of reservoir 
reservoirTop = thickness;  % The top of reservoir 

%% Generate 2D PEBI grid
% 2D PEBI grid is constructed for the provided domain outline, with
% refinement around the wells
nList   = [25]; % Approximate number of cells in x-direction

% geo is a struct that gives the spatial locations of different fractures
geo.xMax = 2250;
geo.yMax = 2250;
geo.xMin = -2250;
geo.yMin = -2250;

geo.fracDepth = [150, 75, -75, -150];

outline = [-2250.0, -2250.0; -2250.0, 2250.0; 2250.0, 2250.0; 2250.0,-2250.0]; % x-y plane size
faceLines = {[0.0,geo.yMax; 0.0,geo.yMin], ...
              [0.0, 150; geo.xMax, 150], ... % horizontal line 1
              [0.0, -75; geo.xMax, -75], ... % horizontal line 2
              [geo.xMin, 75; 0, 75], ... % horizontal line 3
              [geo.xMin, -150; 0, -150] };  % horizontal line 4
 
containFracture = true;
isSplited = true;

if size(faceLines, 1) == 0
    containFracture = false;
end

for i = 1:size(nList, 2)

    %% Make 2.5D reservoir model
    if containFracture
        n   = nList(i); % Approximate number of cells in x-direction
        domainSize   = [ max(outline(:,1)) - min(outline(:,1)), ...
            max(outline(:,2)) - min(outline(:,2)) ];

        [G2D, Pts ] = pebiGrid2D( ...
            max(domainSize)/n, ...
            domainSize          , ...
            'polybdr'        , outline , ...
            'CCRefinement'   , false         , ... % if we need to refine around cells, cellConstraints needs to be given
            'CCFactor'       , 0.1          , ...  % Relative size of well cells
            'CCEps'          , 0.02*max(domainSize) , ...
            'faceConstraints', faceLines    , ...
            'FCRefinement'   , true         , ...  % if we need to refine around faults, faceConstraints needs to be given
            'FCFactor'      , 0.1           , ... % Relative size of fault cells
            'FCEps'         , 0.08*max(domainSize)); % the smaller, the more adrupt the change is

        G2D = removeShortEdges(G2D, 1); % The grid may contain very short edges.
                            % We remove these
    else
        Width = max(outline(:,1));
        Height = max(outline(:,2));
        gS  = min(Width,Height)/nList(i);
        [G2D, Pts ] = pebiGrid2D(gS, [Width,Height]);
    end
    
    if verboseLevel > 1
        fprintf('Total cell number: %d\n', G2D.cells.num);
    end

    wellNo2D                = nan(G2D.cells.num,1);
    wellNo2D(G2D.cells.tag) = 1:nnz(G2D.cells.tag);
    
    G           = makeLayeredGrid(G2D, layerThickness);
    G           = computeGeometry(G);
    G.cells.tag = repmat(G2D.cells.tag, nLayers, 1); % return a repeated copies of array
    %G.faces.tag = false(G.faces.num,1);
    wellNo      = repmat(wellNo2D, nLayers, 1);

    %% Populate with petrophysical properties
    N = [nX , nY, nLayers]; % Num of points in each axial direcion
    K = ones(N) * 1.0; % a pesudo value
    perm = sampleFromBox(G, K, 1:G.cells.num ) * milli * darcy; % CO2 storage  cannot have a very impermeable reservoir
    poro = ones(size(perm))*poroConst;
    rock = makeRock(G, perm, poro);
    
    % The following is to properly selected all cells along fault lines
    faultIndices = [(nLayersInUnderburden+1):(nLayersInUnderburden+nLayersInReservoir)]; % go through the reservoir
    G.cells.tag = repmat(G2D.cells.tag, nLayers, 1); % reservoir layer
    G.faces.tag = repmat(G2D.faces.tag, nLayers, 1);
    [G, G2D] = set_fault_depths(faultIndices, G, G2D, nLayers);
    numFacesAlongFault = length(G.faces.tag(G.faces.tag == 1));

    perm = kron(rock.perm',[1. 1. 1. ]');
    cell_rock_prop = cat(1, perm, rock.poro' , G.cells.volumes');
    
    % save rock properities
    if ~isSplited
        % add pesudo-properties to fracture elements if explicit method is
        % used
        rock.perm = cat(1, rock.perm, zeros(numFacesAlongFault, 1));
        rock.poro = cat(1, rock.poro, zeros(numFacesAlongFault, 1));
    end

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
    
    file_name = strcat( '_res', num2str(n), '_', baseLine);
    vtk_name = strcat( num2str(i, '%0.4d'), file_name );
    fprintf('Writing %s\n', vtk_name);
    vtk_path = fullfile(currentFolder, mesh_folder, vtk_name);
    nPrisms = writeVTKFieldData( G, rock, numFacesAlongFault, ...
                            vtk_path, geo, isSplited);
    
    updateCellBlocksInXml(xmlFilename, nPrisms);
    % output rock properties
    %rock_name = strcat(  num2str(i, '%0.4d'), '_rock.txt' );
    %rock_path = fullfile(currentFolder, rock_folder, rock_name);
    %write_rock_prop(cell_rock_prop, rock_path);

end