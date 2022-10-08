clc
clear all
close all
%given transfer functions
tf1 = tf( 16 , [ 1 3 16 ] );
tf2 = tf( 0.04 , [ 1 0.02 0.04 ] );
tf3 = tf( 1.05*10^7 , [ 1  1600 1.05*10^7 ] );
%experimental results
stepinfo(tf1)
stepinfo(tf2)
stepinfo(tf3)
%results shown on the graphs
figure
step(tf1)
figure
step(tf2)
figure
step(tf3)

