function [G, G2D] = set_fault_depths(faultIndices, G, G2D, nLayers)
    % Loop through each layer
    for i = 1:nLayers
        % Select the range for the current layer
        faceLayerStart = (i - 1) * numel(G2D.faces.tag) + 1;
        faceLayerEnd = i * numel(G2D.faces.tag);

        if ~ismember(i, faultIndices)
            G.faces.tag(faceLayerStart:faceLayerEnd) = false;
        end
    end
end
