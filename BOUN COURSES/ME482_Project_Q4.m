clc
clear
close all

% These will be used in interpolation
beta = [ 1 1.006 1.0246 1.0577 1.1094 1.1867 1.3033 1.4882 1.816 2.5776 ] ;
ratio = 0:0.1:0.9 ;

a = 0.065 ;   % The crack lenght can be adjusted for different problems
K_Ic = 50 ;  % Fracture thoughness given in the problem
d_N = 0.1;  % For better numerical integration

% I will use the Paris-Erdoğan Law to numerically integrate the crack
% lenght. Whenever the calculated K is greater than the fracture
% thoughness, it means that failure has occured.

for N = 1:d_N:1e6
    % Tabular data for the K constant is only upto ratio of 0.9, beyond
    % that interpolation function does not give rea<sonable answers.
    if a >= 0.125*0.9
        disp('a TOO BIG')
        break
    end

    % First, calculate K and using Paris-Erdoğan Law find derivative of a.
    % Then, using Euler integration, find a++.
    K = 40*sqrt(pi*a)*interp1(ratio,beta,a/0.125) ;
    da_dN = 1.1e-11*K^3 ;
    a = a + da_dN*d_N ;

    % This is the failure criterion
    if K >= K_Ic
        disp('SOLUTION CONVERGED')
        a_critical = a ;
        N_fatigue = N ;
        break
    end
end

a_critical
N_fatigue
