function writeVTKHexMesh(G, filename)
    % Get the number of cells and nodes in the grid
    numCells = G.cells.num;
    numNodes = G.nodes.num;
    
    % Create the VTK file
    fileID = fopen(filename, 'w');
    
    % Write the VTK header
    fprintf(fileID, '# vtk DataFile Version 3.0\n');
    fprintf(fileID, 'Hexahedral Mesh\n');
    fprintf(fileID, 'ASCII\n');
    fprintf(fileID, 'DATASET UNSTRUCTURED_GRID\n');
    
    % Write the node coordinates
    fprintf(fileID, 'POINTS %d float\n', numNodes);
    fprintf(fileID, '%f %f %f\n', G.nodes.coords');
    
    % Write the cell connectivity
    fprintf(fileID, 'CELLS %d %d\n', numCells, 9*numCells);
    [nx, ny, nz] = deal(G.cartDims(1), G.cartDims(2), G.cartDims(3));
    for k = 1:nz
        for j = 1:ny
            for i = 1:nx
                %cellIndex = (k-1)*nx*ny + (j-1)*nx + i;
                n1 = (k-1)*(nx+1)*(ny+1) + (j-1)*(nx+1) + i;
                n2 = n1 + 1;
                n3 = n1 + nx + 2;
                n4 = n1 + nx + 1;
                n5 = n1 + (nx+1)*(ny+1);
                n6 = n5 + 1;
                n7 = n5 + nx + 2;
                n8 = n5 + nx + 1;
                fprintf(fileID, '8 %d %d %d %d %d %d %d %d\n', ...
                    n1-1, n2-1, n3-1, n4-1, n5-1, n6-1, n7-1, n8-1);
            end
        end
    end

    
    % Write the cell types (VTK_HEXAHEDRON)
    fprintf(fileID, 'CELL_TYPES %d\n', numCells);
    fprintf(fileID, '%d\n', repmat(12, numCells, 1));
    
    % Write the field data
    fprintf(fileID, 'CELL_DATA %d\n', numCells);
    fprintf(fileID,'SCALARS CELL_MARKERS int 1\n');
    fprintf(fileID,'LOOKUP_TABLE default\n');

    cellCentroids = G.cells.indexMap;
    fieldData = (cellCentroids(:) > 9.0); % using i,j,k logic
    fprintf(fileID, '%d\n', fieldData);
    
    % Close the VTK file
    fclose(fileID);
end