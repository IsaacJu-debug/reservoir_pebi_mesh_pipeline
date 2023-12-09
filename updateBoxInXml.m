function updateBoxInXml(filename, newXMin, newXMax)
    % Read the XML file
    xmlDoc = xmlread(filename);

    % Get all Box elements
    boxElements = xmlDoc.getElementsByTagName('Box');

    % Loop through each Box element
    for i = 0:boxElements.getLength()-1
        box = boxElements.item(i);
        % Check if the name attribute is 'source'
        if strcmp(box.getAttribute('name'), 'source')
            % Update xMin and xMax attributes
            box.setAttribute('xMin', newXMin);
            box.setAttribute('xMax', newXMax);
            break; % Assuming there's only one Box with name='source'
        end
    end

    % Write the modified XML to a new file
    xmlwrite(strrep(filename, '.xml', '_modified.xml'), xmlDoc);
    fprintf('XML file updated and saved as %s\n', strrep(filename, '.xml', '_modified.xml'));
end
