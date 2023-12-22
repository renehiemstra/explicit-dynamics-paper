% initialize
clear all; close all; clc;
addpath('src');

% available methods
Methods = {'standard','fast','lumped'};

%% define spline space
remove_outliers = true;
ndofs  = 100;    % number of elements
bcs = 'fixed';  % (fixed or free) boundary conditions 

colors = [ [0.4660, 0.6740, 0.1880];
        [0, 0.4470, 0.7410]  ;
        [0.8500, 0.3250, 0.0980];
        [0.9290, 0.6940, 0.1250];
        [0.4940, 0.1840, 0.5560];	  
        [0.3010, 0.7450, 0.9330];
        [0.6350, 0.0780, 0.1840]];

p = 4;

[Rw1, Rx1] = spectrum_and_modes_bar('standard', p, ndofs, bcs, remove_outliers==true);
[Rw2, Rx2] = spectrum_and_modes_bar('fast', p, ndofs, bcs, remove_outliers==true);
[Rw3, Rx3] = spectrum_and_modes_bar('lumped', p, ndofs, bcs, remove_outliers==true);


%% Plot results

fontsize = 18;
m1 = length(Rw1);
m2 = length(Rw2);

% create figures and axes
fig1 = figure('Position',[100,100,700,450]); set(gcf, 'Color', 'w');
fig2 = figure('Position',[100,100,700,450]); set(gcf, 'Color', 'w');
axis1 = axes('Parent',fig1); hold(axis1,'on');
axis2 = axes('Parent',fig2); hold(axis2,'on'); 

% Normalized eigen frequencies
plot((1:m2)./m2, ones(m2,1), 'k', 'parent',axis1,'linewidth',0.25, 'DisplayName', 'Exact');
plot((1:m1)./m2, Rw1,'-','parent',axis1,'Color','k', 'linewidth', 2, 'DisplayName', 'Consistent mass');
plot((1:m2)./m2, Rw2,'-','parent',axis1,'Color', colors(2,:), 'linewidth', 1, 'DisplayName', 'Approximate mass');
plot((1:m2)./m2, Rw3,'-','parent',axis1,'Color', colors(3,:), 'linewidth', 1, 'DisplayName', 'Lumped mass');

% figure and axis properties
set(axis1, 'ylim', [0,1.1]);
set(axis1,'XScale', 'linear', 'xlim', [0, 1], 'YScale', 'linear', 'ylim', [0, 1.1]);
xlabel(axis1, '$n / N$','interpreter','latex', 'FontSize',fontsize);
ylabel(axis1, '$\omega^h / \omega$','interpreter','latex', 'FontSize',fontsize);
legend1 = legend(axis1,'show');
set(legend1,'Position',[0.687857142857143 0.245555555555556 0.167142857142857 0.0844444444444443]);

% semilogy-errors in frequencies
plot((1:m1)./m2, abs(Rw1-1),'-','parent',axis2,'Color','k', 'linewidth', 2, 'DisplayName', 'Consistent mass');
plot((1:m1)./m2, abs(Rw2-1),'-','parent',axis2,'Color',colors(2,:), 'linewidth', 1, 'DisplayName', 'Approximate mass');
plot((1:m1)./m2, abs(Rw3-1),'-','parent',axis2,'Color',colors(3,:), 'linewidth', 1, 'DisplayName', 'Lumped mass');

% figure and axis properties
set(axis2,'XScale', 'linear', 'xlim', [0, 1], 'XTick', 0.0:0.2:1.0, 'YScale', 'log', 'ylim', [1e-15, 1.1], 'YTick', [1e-15, 1e-12, 1e-9, 1e-6, 1e-3, 1]);
xlabel(axis2, '$n / N$','interpreter','latex', 'FontSize',fontsize);
ylabel(axis2, '$| \omega^h - \omega | / \omega$','interpreter','latex', 'FontSize',fontsize);
legend2 = legend(axis2,'show');
set(legend2,'Position',[0.687857142857143 0.245555555555556 0.167142857142857 0.0844444444444443]);