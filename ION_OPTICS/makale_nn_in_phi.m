clc
clear 
close all

load custom3.mat

N_radial = 151 ;
N_axial  = 376 ;
N_iteration = 250;

V_print_type = "contour";
print_grid = false;
grid_color = [0.9290 0.6940 0.1250];

% camera_location = "TOP" ;
% camera_location = "LEFT" ;
camera_location = "FRONT" ;
% camera_location = "ISOMETRIC" ;

name1 = "D:\C++ PROJECTS\RESULTS\001 (375)\[res]" ;

smooth_ratio = 5;
print_step = 5 ;
plot_step = 1 ;

% type = "phi" ;  
type = "in-time" ; 
% type = "en-time" ;
% type = "nn-time" ;
% type = "RHO-time" ;
% type = "Qin-time" ;

if N_radial ==101
    x1 = [51 71 71 51] ;
    y1 = [51 51 101 101] ;
    x2 = [130 171 171 130] ;
    y2 = [30 30 101 101] ;
elseif N_radial == 51
    x1 = [51 71 71 51] ;
    y1 = [51 51 101 101] ;
    x2 = [130 171 171 130] ;
    y2 = [30 30 101 101] ;
elseif N_radial == 151
    x1 = [51 71 71 51] ;
    y1 = [51 51 101 101] ;
    x2 = [130 171 171 130] ;
    y2 = [30 30 101 101] ;
end

z = [10^14 10^14 10^14 10^14] ;

if type == "phi"
    type_sol = "ELECTRIC POTENTIAL" ;
elseif type == "in-time"
    type_sol = "ION NUMBER DENSITY" ;
elseif type == "nn-time"
    type_sol = "NEUTRAL NUMBER DENSITY" ;
elseif type == "en-time"
    type_sol = "ELECTRON CHARGE DENSITY";
elseif type == "RHO-time"
    type_sol = "CHARGE DENSITY";
elseif type == "Qin-time"
    type_sol = "ION CHARGE DENSITY";
end

name2 = "-iteration" ;

for i = print_step:print_step:(N_iteration-1)
    number = num2str(i) ;
    file_name = sprintf("%s%s%s%s.txt",name1,type,name2,number) ;
    V(i/print_step,:) = get_data_tex(file_name, [1 N_axial*N_radial] ) ;
end

for i = 1:(N_iteration-1)/print_step
    for radial = 1:N_radial
        for axial = 1:N_axial
            OneDto2D(i,radial,axial) = V(i,(axial-1)*N_radial + radial) ;
        end
    end
end

f=figure;
f.WindowState="fullscreen";

for i = 1:(N_iteration-1)/print_step  
    for k = 1:N_radial
        for j = 1:N_axial
            if (k <= (smooth_ratio-1)/2) || (k >= (N_radial - (smooth_ratio-1)/2)) || (j <= (smooth_ratio-1)/2) || (j >= (N_axial - (smooth_ratio-1)/2))
                One(k,j) = OneDto2D(i,k,j);
            else
                SUM = sum( sum( OneDto2D( i, (k-(smooth_ratio-1)/2):(k+(smooth_ratio-1)/2), (j-(smooth_ratio-1)/2):(j+(smooth_ratio-1)/2 )) ) );
                One(k,j) = SUM/smooth_ratio^2;
            end
        end
    end

    hold off

    if type == "phi" && V_print_type == "contour"
        contour(1:plot_step:N_axial, 0:plot_step:N_radial-1, One,"fill",1,"LevelStep",100)
        hold on
        contour(1:plot_step:N_axial, 0:-plot_step:-N_radial+1, One,"fill",1,"LevelStep",100)
%         colormap("colorcube")
        colorbar
    else
        surf(1:plot_step:N_axial, 0:plot_step:N_radial-1, One)
        hold on
        surf(1:plot_step:N_axial, 0:-plot_step:-N_radial+1, One)
        colormap(custom3)
%         colormap("colorcube")
        colorbar
    end
    if print_grid
        hold on
        fill3(x1, y1, z,grid_color)
        hold on
        fill3(x2, y2, z,grid_color)
        hold on
        fill3(x1, -y1, z,grid_color)
        hold on
        fill3(x2, -y2, z,grid_color)
    end
 
    title_name = sprintf("%s - TIME STEP : %d", type_sol, i*print_step) ;
    title(title_name)
%     axis equal
    xlim([0 N_axial])
    ylim([-N_radial N_radial])

    if camera_location == "FRONT"
        view([500 0 0])   
    elseif camera_location == "TOP"
        view([0 0 500])    
    elseif camera_location == "LEFT"
        view([0 -500 500])  
    elseif camera_location == "ISOMETRIC"
        view([500 -500 500]) 
    end
    
%     pause(0.5)

    drawnow
end


