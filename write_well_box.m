function [wellCellIDs] = write_well_box( iCase, G, wellNo, filename, ...
                                    reservoirBottom, reservoirTop )
% filename: basefile where source information is given

well_markers = zeros(G.cells.num,1);
well_markers( wellNo > 0 ) = 1;
cellIDs = 1:G.cells.num;
wellCellIDs = cellIDs(well_markers > 0);
if length(wellCellIDs) ~= 1
  error("Only one well cell at the moment");
end

numFacesPerCell = diff(G.cells.facePos);
nFaces = numFacesPerCell(wellCellIDs);

cellFaces = G.cells.faces(G.cells.facePos(wellCellIDs) : G.cells.facePos(wellCellIDs+1)-1, 1);

tmp = [];
for iFace = 1:nFaces
  faceID = cellFaces(iFace);
  nodes = G.faces.nodes(G.faces.nodePos(faceID) : G.faces.nodePos(faceID+1)-1, :) - 1;
  tmp = [tmp; nodes];
end
tmp = unique(tmp)+1;

xMax = max(G.nodes.coords(tmp,1));
xMin = min(G.nodes.coords(tmp,1));
yMax = max(G.nodes.coords(tmp,2));
yMin = min(G.nodes.coords(tmp,2));
zMax = max(G.nodes.coords(tmp,3));
zMin = min(G.nodes.coords(tmp,3));

hx = xMax-xMin;
hy = yMax-yMin;
hz = reservoirTop-reservoirBottom;
tolX = 1.e-6*hx;
tolY = 1.e-6*hy;
tolZ = 1.e-6*hz;
FID_box = fopen(strcat( num2str(iCase,'%0.4d'),'_box_source.xml'),'w');

fprintf(FID_box,'%s\n','<?xml version="1.0" ?>');
fprintf(FID_box,'\n');
fprintf(FID_box,'%s\n','<Problem>');
fprintf(FID_box,'%s\n','  <Geometry>');
fprintf(FID_box,'%s\n','    <Box');
fprintf(FID_box,'%s\n','      name="source"');
fprintf(FID_box,'%s'  ,'      xMin="{ ');
fprintf(FID_box,'%13.7e', xMin-tolX);
fprintf(FID_box,'%s'  ,', ');
fprintf(FID_box,'%13.7e', yMin-tolY);
fprintf(FID_box,'%s'  ,', ');
fprintf(FID_box,'%13.7e', zMin-tolZ);
fprintf(FID_box,'%s\n'  ,' }"');
fprintf(FID_box,'%s'  ,'      xMax="{ ');
fprintf(FID_box,'%13.7e', xMax+tolX);
fprintf(FID_box,'%s'  ,', ');
fprintf(FID_box,'%13.7e', yMax+tolY);
fprintf(FID_box,'%s'  ,', ');
fprintf(FID_box,'%13.7e', zMax+tolZ);
fprintf(FID_box,'%s\n'  ,' }"/>');
fprintf(FID_box,'%s\n','  </Geometry>');
fprintf(FID_box,'%s\n','</Problem>');

% Calculate new xMin and xMax
newXMin = sprintf('{ %f, %f, %f }', xMin - tolX, yMin - tolY, reservoirBottom - tolZ);
newXMax = sprintf('{ %f, %f, %f }', xMax + tolX, yMax + tolY, reservoirTop + tolZ);

fclose(FID_box);

% Read the XML file
updateBoxInXml(filename, newXMin, newXMax);


end

