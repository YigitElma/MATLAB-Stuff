clc
clear all
close all

A=[5 -1 1; 2 8 -1; -1 1 4];
B=[10; 11; 3];

P=[0; 0; 0];

[X_jacobi k_jacobi]=jacobi(A, B, P)