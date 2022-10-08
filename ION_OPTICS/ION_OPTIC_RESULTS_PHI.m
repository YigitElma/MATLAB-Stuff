clc
clear
close all

load imported_data.mat

for i = 1:10:909
    One(:,:) = OneDto2D(i,1:5:700,:) ;

    surf(1:50, 1:5:700, One)
    colormap("hot")
    colorbar
    title("iteration number #", i*10 )
    view([0 0 500])

    drawnow
end






