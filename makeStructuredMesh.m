clear
close
%% Define the dimensions of the grid
nx = 3; % Number of cells in the x-direction
ny = 3;  % Number of cells in the y-direction
nz = 2;  % Number of cells in the z-direction

% Define the physical dimensions of the grid
xmin = 0; xmax = 1; % Range for x
ymin = 0; ymax = 1; % Range for y
zmin = 0; zmax = 2; % Range for z

% Create the structured grid
G = cartGrid([nx, ny, nz], [xmax-xmin, ymax-ymin, zmax-zmin]);
G.nodes.coords = G.nodes.coords + repmat([xmin, ymin, zmin], G.nodes.num, 1);

%% Visualize the grid
figure;
plotGrid(G, 'FaceColor', 'r', 'EdgeColor', 'k', 'FaceAlpha', 0.5);
view(3);
axis equal tight;
xlabel('x');
ylabel('y');
zlabel('z');
title('Customized Structured Grid');

%% Write VTK files
folder = '01_mesh';
filename = 'hexa_two_body.vtk';

% Create the folder if it doesn't exist
if ~exist(folder, 'dir')
    mkdir(folder);
end

% Generate the VTK file
fullPath = fullfile(folder, filename);
write_hexa_vtk(G, fullPath);
