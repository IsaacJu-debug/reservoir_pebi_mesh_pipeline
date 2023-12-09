function write_vtk( G, Pts, rock, attribute, wellNo, filename )

  %% Allocating space for polyhedra output
  numFacesPerCell = diff(G.cells.facePos);
  numNodesPerFace = diff(G.faces.nodePos);
  
  nPrisms = zeros(9,1 ); %prism3 (aka wedge), ..., prism11
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
    
  % % Well cell bounding box
  well_markers = zeros(G.cells.num,1);
  well_markers( wellNo > 0 ) = 1;

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
  
  permeability = kron(rock.perm',[1. 1. 1. ]');
  permeability(3,:) = 0.1*permeability(3,:);
  fprintf(FID_vtk,'VECTORS PERM float\n');
  fprintf(FID_vtk, ...
          '%11.5e %11.5e %11.5e\n', ...        
          permeability );     

%  fprintf(FID_vtk,'VECTORS PRESSURE_POINTS float\n');
%  fprintf(FID_vtk, ...
%          '%11.5e\n', ...        
%          [Pts, G.cells.centroids(:,3)]' );           
        
  fprintf(FID_vtk,'SCALARS PORO float 1\n');
  fprintf(FID_vtk,'LOOKUP_TABLE default\n');
  fprintf(FID_vtk, ...
          '%11.5e\n', ...
          rock.poro); 

  fprintf(FID_vtk,'SCALARS REGION float 1\n');
  fprintf(FID_vtk,'LOOKUP_TABLE default\n');
  fprintf(FID_vtk, ...
          '%d\n', ...
          attribute); 
  
%  fprintf(FID_vtk,'SCALARS WELL_MARKERS int 1\n');
%  fprintf(FID_vtk,'LOOKUP_TABLE default\n');
%  fprintf(FID_vtk, ...
%          '%d\n', ...
%          well_markers );       
        
  fclose( FID_vtk);
  
end 
