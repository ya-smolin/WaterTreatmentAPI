clear all; close all; clc;
R_c = @(k, C)C.*(1+2*sqrt(k(1)/k(2)))./(k(1)+C+C.^2/k(2));%looped

N = 100000;
k1 = linspace(1, 100, N);
k2 = linspace(1, 100, N);
C = 1:1:10;

tic;
for i = 1:N
    R_cal([k1(i) k2(i)], C);
end
toc;

tic;
for i = 1:N
    R_c([k1(i) k2(i)], C);
end
toc;

