clear all
close all
clc
%%

syms a1 a2 psi z

C=[0 a1 0 0 0; a1 0 a2 0 0; 0 a2 0 a2 0; 0 0 a2 0 a1; 0 0 0 a1 0];

psi0=[1,0,0,0,0]';
psi=expm(-i*C*z)*psi0;

%%
syms psii t

v=2;
H=[0 v 0;v 0 v;0 v 0];
psi00=[0,1,0]';
psii=expm(-i*H*t)*psi00