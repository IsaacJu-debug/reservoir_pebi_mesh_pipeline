function nPrisms = write_vtk_field_data( G, rock, numFacesAlongFault, filename, ...
                                        reservoirBottom, reservoirTop, isSplited)
  
  %% allocating space for faulted faces
  %numFacesAlongFault = length(G.faces.tag(G.faces.tag == 1));
  if isSplited
    numFacesToBeSplited = numFacesAlongFault; % only split faces
    numFacesAlongFault = 0;
  else
    numFacesToBeSplited =0; % no fault markers for splitting faces
  end

  %% Allocating space for polyhedra output
  numFacesPerCell = diff(G.cells.facePos);
  numNodesPerFace = diff(G.faces.nodePos);
  
  nPrisms = zeros(9,1 ); %prism3 (aka wedge), prism4 (hexhedra), ..., prism11
  for nSides = 3:11
    nPrisms(nSides-2) = sum(numFacesPerCell == ( nSides + 2) );
  end
  
  % GEOSX supports prisms with at most 11-edge sided polygonal basois
  if sum(nPrisms) ~= G.cells.num
    error("Cells are not just prism3 to prism11");
  end
  
  storage = 0;
  for nSides = 3:11
    nFaces = nSides + 2;
    offset = nFaces + 1 + 2*nSides + (nFaces-2)*4;
    storage = storage + nPrisms(nSides-2)*offset;
  end
  
  %% Populate connectivity array
  % Connectivity array contains two part: polygons and quad. All quad cells
  % are faulted faces
  faultFaceEntries = numFacesAlongFault * 5;
  connectivity = zeros(storage+G.cells.num + faultFaceEntries,1);
  
  ptr = 1;
  for iCell = 1:G.cells.num
    
    nFaces = numFacesPerCell(iCell);
    nSides = nFaces - 2;
    offset = nFaces + 1 + 2*nSides + (nFaces-2)*4;
    
    connectivity( ptr ) = numFacesPerCell(iCell);
    
    tmp = [ offset; nFaces ];
    % cells.faces stores cell_to_face map
    % cells.facePos stores index of faces belonging to the current cell
    cellFaces = G.cells.faces(G.cells.facePos(iCell) : G.cells.facePos(iCell+1)-1, 1);
    
    for iFace = 1:nFaces
      faceID = cellFaces(iFace);
      N = numNodesPerFace( faceID );
      % faces.nodes stores face_to_vertex map
      nodes = G.faces.nodes(G.faces.nodePos(faceID) : G.faces.nodePos(faceID+1)-1, :) - 1;
      tmp = [tmp; N; nodes];
  
    end
    
    connectivity(ptr:ptr+offset) = tmp;
    
    ptr = ptr + offset+1;  
  
  end
  %% cell region attributes
  cell_markers = ones(G.cells.num + numFacesAlongFault,1); % reservoir
  cell_markers( G.cells.centroids(:,3) < reservoirBottom ) = 2; % underburden
  cell_markers( G.cells.centroids(:,3) > reservoirTop ) = 3; % overburden
  if ~isSplited
    cell_markers(G.cells.num+1:G.cells.num + numFacesAlongFault) = 4; % Fault faces
  end

  %% fault faces
  % First, locate all faulted faces. Treat these faces as a quad cell,
  % appended to cell block. CELL_TYPES need to be updated; CELL_DATA fields
  % need to be updated. 
  % face.tages gives fault faces
  % use face.nodes (face_to_vertex) to locate corresponding vertices 
  
  numExtraFace = length(G.faces.neighbors) - length(G.faces.tag);
  faceArray = find(G.faces.tag == 1) + numExtraFace; % 1614 corresponds to extra faces added because of mesh extruding.                                                                    
  faceOffset = 4;

  for iFace = 1:numFacesAlongFault
      
    tmp = [faceOffset];
    faceID = faceArray(iFace); 
    % faces.nodes stores face_to_vertex map
    nodes = G.faces.nodes(G.faces.nodePos(faceID) : G.faces.nodePos(faceID+1)-1, :) - 1;
    %tmp = [tmp; nodes];
    tmp = [tmp; nodes];
    connectivity(ptr:ptr+faceOffset) = tmp;
    ptr = ptr + faceOffset + 1;
    % assign cells along fault to be faulted cells

  end

  if isSplited
    % adding extra field for identifying fractured cells
    fault_markers = zeros(G.cells.num, 1); % rest of domain 
    for iFace = 1:numFacesToBeSplited
    
        faceID = faceArray(iFace); 
        fault_markers(G.faces.neighbors(faceID, 1)) = 100; % one side
        fault_markers(G.faces.neighbors(faceID, 2)) = 101; % another side
        
    end
  end
  
  cellTypeArray = 9*ones(G.cells.num + numFacesAlongFault,1); % all quad
  cellTypeArray(1:G.cells.num) = 42; % change to 42
  %% OUTPUT vtk
  % Open output file
  FID_vtk = fopen(filename,'w');
  
  % Header
  fprintf(FID_vtk,'%s\n','# vtk DataFile Version 2.0');
  fprintf(FID_vtk,'%s\n','pebi grid mesh');
  fprintf(FID_vtk,'%s\n','ASCII');
  fprintf(FID_vtk,'%s\n','DATASET UNSTRUCTURED_GRID');
  
  fprintf(FID_vtk,'POINTS %d float\n',G.nodes.num);
  fprintf(FID_vtk, ...
          '%13.7e %13.7e %13.7e\n', ...
          G.nodes.coords' );
         
  fprintf(FID_vtk,'CELLS %d %d\n',G.cells.num + numFacesAlongFault, length(connectivity)  );
  fprintf(FID_vtk, ...
          '%d %d %d %d %d %d %d %d %d %d\n', ...
          connectivity );      
  fprintf(FID_vtk, '\n');      
        
  fprintf(FID_vtk,'CELL_TYPES %d\n',G.cells.num + numFacesAlongFault);
  fprintf(FID_vtk, ...
          '%d\n', ...
          cellTypeArray );      
  fprintf(FID_vtk, '\n');      
    
  %% cell data
  fprintf(FID_vtk,'CELL_DATA %d\n', length(cell_markers));  

  fprintf(FID_vtk,'SCALARS CELL_MARKERS int 1\n');
  fprintf(FID_vtk,'LOOKUP_TABLE default\n');
  fprintf(FID_vtk, ...
          '%d\n', ...
          cell_markers );  

  %% field data
  cell_size = length(cell_markers);
  field_number = 2;
  if isSplited
      field_number = 3; % we need to add one more field for hosting fault
  end
  fprintf(FID_vtk, '\n');  
  fprintf(FID_vtk,'FIELD FieldData %d\n', field_number);  % the number should correspond to the number of fields
  fprintf(FID_vtk,'PERM 3 %d float\n', cell_size);
  fprintf(FID_vtk, ...
          '%11.5e %11.5e %11.5e\n', ...        
          kron(rock.perm',[1. 1. 1. ]') );     
       
%   presPointsArray = cat(2, [Pts, G.cells.centroids(:,3)]', zeros(3,numFacesAlongFault));
%   fprintf(FID_vtk,'PRESSURE_POINTS 3 %d float\n', cell_size);
%   fprintf(FID_vtk, ...
%           '%11.5e %11.5e %11.5e\n', ...        
%           presPointsArray );           
        
  fprintf(FID_vtk,'PORO 1 %d float\n', cell_size);
  %fprintf(FID_vtk,'LOOKUP_TABLE default\n');
  fprintf(FID_vtk, ...
          '%11.5e\n', ...
          rock.poro); 
  
  if isSplited
    %% cell data
    fprintf(FID_vtk,'FAULT 1 %d int\n', cell_size);
    fprintf(FID_vtk, ...
          '%d\n', ...
          fault_markers ); 
  end
  fclose( FID_vtk);
  
end 
