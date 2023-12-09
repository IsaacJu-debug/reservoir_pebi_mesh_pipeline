function savePcolor(kMap, dir, name)
%SAVEPCOLOR Summary of this function goes here
%   Detailed explanation goes here
    if ~exist(dir, 'dir')
        mkdir(dir);
        disp(strcat('mkdir ', dir, '...'));
    end
    f = figure('visible','off');
    colormap("jet");
    pcolor(kMap);
    colorbar;
    ax = gca;
    path = strcat(dir, '/', name);
    exportgraphics(ax, path);

end

