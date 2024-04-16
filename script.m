%% initialize
clear all; close all; clc;
addpath('src');

% spline discretization
ne = 1;
p = 3;
a = 0; b = 1.0;
% open knot vector
kts  = [a*ones(1,p) linspace(a, b, ne+1) b*ones(1,p)];
n = length(kts) - p - 1;

% contol points of the mapping
alpha = 1;
gp = grevillepts(p, kts);
[xcpts, ycpts] = meshgrid(gp,gp); 
xcpts(2,2:3) = xcpts(2,2:3) + alpha;
xcpts(3,2:3) = xcpts(3,2:3) - alpha;
ycpts(2:3,2) = ycpts(2:3,2) + alpha;
ycpts(2:3,3) = ycpts(2:3,3) - alpha;
z = zeros(n,n);

% compute mapping
mm = 20;
uu = linspace(0.0,1.0,mm);
[X, dX, ddX] = mapping(p, kts, uu, xcpts, ycpts);

% plot mesh
mesh(X{1}, X{2}, zeros(mm,mm));