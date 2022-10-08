% -------------------------------------------------------------------------
% 
% Yiğit Günsür ELMACIOĞLU  -  2017405120
% 
% -------------------------------------------------------------------------

clc
clear 
close all

k = 1.4;
p_exit_req = 0.7 ;
tolerance = 0.0001 ;

Area = @(x) 0.1 + x^2;
P_0_after_shock = @(M) ((k+1)/2*M^2)^(k/(k-1)) * (1+(k-1)/2*M^2)^( -k/(k-1) ) / ( 2*k/(k+1)*M^2 - (k-1)/(k+1) )^(1/(k-1));
M_after_shock = @(M) sqrt( (M^2 + (2/(k-1)) ) / (2*k/(k-1)*M^2 - 1) );

A_throat = Area(0);
[p_in, M_in] = isentropic_ratio_sub(Area(-0.5)/A_throat);

% %  This part is for PART a
% x_shock_04 = 0.4;
% [p_shock_1, M_1] = isentropic_ratio_sup(Area(x_shock_04)/A_throat)
% p0_shock_2 = P_0_after_shock(M_1) 
% M_2 = M_after_shock(M_1)
% A_throat_new = Area(x_shock_04)*M_2 / ( 2/(k+1) * (1+(k-1)/2*M_2^2) )^((k+1)/(2*k-2))
% [p_exit_rel, M_exit] = isentropic_ratio_sub(Area(0.5)/A_throat_new)
% p_exit_04 = p_exit_rel * p0_shock_2 

for x_shock = 0.01:tolerance:0.5

    [p_shock_1, M_1] = isentropic_ratio_sup(Area(x_shock)/A_throat);
    p0_shock_2 = P_0_after_shock(M_1) ;
    M_2 = M_after_shock(M_1);
    A_throat_new = Area(x_shock)*M_2 / ( 2/(k+1) * (1+(k-1)/2*M_2^2) )^((k+1)/(2*k-2));
    [p_exit_rel, M_exit] = isentropic_ratio_sub(Area(0.5)/A_throat_new);
    p_exit = p_exit_rel * p0_shock_2 ;

    if abs(p_exit - p_exit_req) < tolerance
        x_shock_normal = x_shock ;
        break
    end

end

x_shock_normal

% there are several roots of the following equation. Because of this, I
% used two functions with different initial value for fzero(), one for
% subsonic speeds and one for supersonic speeds.
function [p0_ratio,M] = isentropic_ratio_sub(A)
k = 1.4;
func = @(M) 1/M * ( 2/(k+1) * (1+(k-1)/2*M^2) )^((k+1)/(2*k-2)) - A;
M = fzero(func,1);              
p0_ratio = (1+(k-1)/2*M^2)^( -k/(k-1) );
end

function [p0_ratio,M] = isentropic_ratio_sup(A)
k = 1.4;
func = @(M) 1/M * ( 2/(k+1) * (1+(k-1)/2*M^2) )^((k+1)/(2*k-2)) - A;
M = fzero(func,2);
p0_ratio = (1+(k-1)/2*M^2)^( -k/(k-1) );
end