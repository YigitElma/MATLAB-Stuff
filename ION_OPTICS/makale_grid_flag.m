clc
clear 
close all

tic

N_radial = 150 ;
N_axial  = 375 ;
N_theta = 10 ;

N_iteration = 51 ;

plot_step = 1 ;

name = "D:\C++ PROJECTS\PIC-DSMC\gridflags.txt" ;
V(:) = get_data_tex(name, [(N_theta-1)*N_axial*N_radial+1 N_theta*N_axial*N_radial] ) ;


for radial = 1:N_radial
    for axial = 1:N_axial
        OneDto2D(radial,axial) = V((radial-1)*N_axial + axial) ;
    end
end


surf(1:plot_step:N_axial, 1:plot_step:N_radial, OneDto2D(1:plot_step:N_radial, 1:plot_step:N_axial))
colorbar
title(" GEOMETRY and CELL TYPES")
view([0 0 500])


toc