close all;
clear all;
clc;

syms k t u
f = (k.^2*t.^2-2*u*t+u.^2)./(2*t);

diff(f,u)