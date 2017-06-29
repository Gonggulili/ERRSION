clc
clear all
close all

Nr=32;
U=8;
H=wgn(Nr,U,0,'complex');
A=H'*H;
L=outCholDecompose(A);
A_=L*L';
A_-A

R=chol(A)