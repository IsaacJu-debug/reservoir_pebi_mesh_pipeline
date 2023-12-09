function write_well_index( iCase, G, rock, wellCellIDs )

W = addWell([], G, rock, wellCellIDs );
wellIndex = W.WI;


FID_well = fopen(strcat( num2str(iCase,'%0.4d'),'_well_info.txt'),'w');



for iWell = 1:length( wellCellIDs )
  ID = wellCellIDs(iWell);
  fprintf( FID_well, '%s\n', '# CELL_ID    X_CENTROID    Y_CENTROID    Z_CENTROID          PERM          PORO    WELL_INDEX');
  fprintf( FID_well, ...
           '% 9d %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e\n', ...
           wellCellIDs( iWell ), ...
           G.cells.centroids( ID, :), ...
           rock.perm( ID), ...
           rock.poro( ID), ...
           wellIndex( iWell ) );
end
         
fclose(FID_well);

end

