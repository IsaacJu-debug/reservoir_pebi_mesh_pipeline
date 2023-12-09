function write_xml( nPrisms, filename, baseName, vtkname)

    cell_block_list = [ "wedges"; "hexahedra"; "pentagonalPrisms"; "hexagonalPrisms";
                        "heptagonalPrisms";  "octagonalPrisms";  "nonagonalPrisms";
                         "decagonalPrisms";  "hendecagonalPrisms"];
    
    baseStruct = readstruct(baseName);
    baseStruct.Mesh.VTKMesh.fileAttribute = convertCharsToStrings(vtkname);
    
    cell_array = cell_block_list(nPrisms > 0);
    for iCell = 0:1
        newLine =  '{';
        for iBlock = 1:length(cell_array)
            tmp = sprintf(' %d_%s', iCell, cell_array(iBlock));
            newLine = append(newLine, tmp); 
            if iBlock < length(cell_array)
                newLine = append(newLine, ',');
            end
        end
        newLine = append(newLine, ' }');
        baseStruct.ElementRegions.CellElementRegion(iCell+1).cellBlocksAttribute = convertCharsToStrings(newLine);
    end
    
    writestruct(baseStruct, filename, "StructNodeName", 'Problem');
