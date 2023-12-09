function updateCellBlocksInXml(filename, nPrisms)
    % Define the list of cell types
    cell_block_list = [ "wedges"; "hexahedra"; "pentagonalPrisms"; "hexagonalPrisms";
                    "heptagonalPrisms";  "octagonalPrisms";  "nonagonalPrisms";
                     "decagonalPrisms";  "hendecagonalPrisms"];
    cell_array = cell_block_list(nPrisms > 0);

    % Read the XML file
    xmlDoc = xmlread(filename);

    % Find all <CellElementRegion> nodes
    cellElementRegions = xmlDoc.getElementsByTagName('CellElementRegion');

    % Loop through each <CellElementRegion> and update cellBlocks attribute
    for i = 0:cellElementRegions.getLength()-1
        region = cellElementRegions.item(i);

        % Construct the cellBlocks string
        tmp = strjoin(arrayfun(@(x, y) sprintf('%d_%s', i+1, y), nPrisms(nPrisms > 0), cell_array, 'UniformOutput', false), ', ');
        % Update cellBlocks attribute
        tmp = ['{', tmp,'}']
        region.setAttribute('cellBlocks', tmp);
    end

    % Write the modified XML to a new file
    newFilename = strrep(filename, '.xml', '_modified.xml');
    xmlwrite(newFilename, xmlDoc);
    fprintf('XML file updated and saved as %s\n', newFilename);
end
