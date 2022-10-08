%ME 335 Homework 1 Question 4
%02.05.2021

clc
clear all
close all

K = 11.43;      %can be calculated by the formula
tau = 0.081;    %can be calculated by the formula

sys = tf( [K], [tau 1] );   %constructing transfer function

h = 0.001;    %this is to increase sharpness of the step response
t1 = 0:h:0.25-h;   %between this period, input is 1
N = length(t1);
t2 = 0.25:h:0.5;   %between this period, input is 5
M = length(t2);
t3 = 0.5+h:h:0.6;  %between this period, input is 1
L = length(t3);

U = [ ones(1,N) 5*ones(1,M) ones(1,L) ];   %I combined all to have only one input vector
t = 0:h:0.6;

lsim(sys,U,t)  %this function calculates and plots response of transfer function to input vector
grid on
