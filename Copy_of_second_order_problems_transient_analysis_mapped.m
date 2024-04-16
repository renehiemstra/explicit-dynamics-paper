%% initialize
clear all; close all; clc;
addpath('src');

export = false;
plt_on = false;

%% Settings
D = @(r,t,theta) 0.5*analytical_solution_annulus_displacement(r,t,theta);
V = @(r,t,theta) 0.5*analytical_solution_annulus_velocity(r,t,theta);

% define spline space
P   = [3,4, 5];                                     % polynomial degree
N   = [8,16,32];                                    % number of elements
Tmax = 2*pi / 1.106470948850117e+1;

% allocate space for residual

R =     {zeros(length(N),length(P)) zeros(length(N),length(P)) zeros(length(N),length(P)); 
         zeros(length(N),length(P)) zeros(length(N),length(P)) zeros(length(N),length(P))};
T =     {zeros(length(N),length(P)) zeros(length(N),length(P)) zeros(length(N),length(P));
         zeros(length(N),length(P)) zeros(length(N),length(P)) zeros(length(N),length(P))};
Ndofs = {zeros(length(N),length(P)) zeros(length(N),length(P)) zeros(length(N),length(P));
         zeros(length(N),length(P)) zeros(length(N),length(P)) zeros(length(N),length(P))};

% loop over degrees and elements
M = length(P) * length(N);
m = 1; l = 1;
for method={'standard','fast','lumped'}
    k = 1;
    for remove_outliers=[0,1]
        for j=1:length(P)
           for i=1:length(N)
                % compute solution to linear transient problem
                p = P(j);
                n = N(i);
                tfac = -1; % ((p/2)./n).^p;
                [l2e, dt_crit, frame1, frame2, ndofs]  = transient_analysis_annulus(D, V, p, n, Tmax, tfac, method, remove_outliers, plt_on);
    
                % error
                R{k,l}(i,j) = l2e;
                T{k,l}(i,j) = dt_crit;
                Ndofs{k,l}(i,j) = ndofs;
            
                sprintf('Case %d%s%d%s%d%s%d ', m, '/', M, ': degree ', p, ' ndofs ', n)
                m = m+1;
           end
        end
        k = k+1;
    end
    l = l+1;
end

%% plot results
close all;

export = true;

colors = [ [0.4660, 0.6740, 0.1880];
        [0, 0.4470, 0.7410]  ;
        [0.8500, 0.3250, 0.0980];
        [0.9290, 0.6940, 0.1250];
        [0.4940, 0.1840, 0.5560];	  
        [0.3010, 0.7450, 0.9330];
        [0.6350, 0.0780, 0.1840]];

% create figures and axes
fig1  = figure('Position',[0,0,600,400], 'name','Convergence without using outlier removal'); set(gcf, 'Color', 'w');
axis1 = axes('Parent',fig1); hold(axis1,'on');

% create figures and axes
fig2  = figure('Position',[0,0,600,400], 'name','Convergence using outlier removal'); set(gcf, 'Color', 'w');
axis2 = axes('Parent',fig2); hold(axis2,'on');

% Normalized eigen f
for k=1:length(P)
    p = P(k);
    dname = {sprintf('%s%d', 'Consistent mass p=', p), sprintf('%s%d', 'Approximate mass p=', p), sprintf('%s%d', 'Lumped mass p=', p)};
    % plot results without outlier removal
    plot(sqrt(Ndofs{1,1}(:,k)),R{1,1}(:,k),'-','DisplayName',dname{1},'linewidth',2,'Color',colors(k,:),'MarkerSize',10,'MarkerFaceColor', 'w','parent',axis1); 
    plot(sqrt(Ndofs{1,2}(:,k)),R{1,2}(:,k),'+','DisplayName',dname{2},'linewidth',2,'Color',colors(k,:),'MarkerSize',10,'MarkerFaceColor', 'w','parent',axis1); 
    plot(sqrt(Ndofs{1,3}(:,k)),R{1,3}(:,k),'--','DisplayName',dname{3}, 'linewidth',1,'Color',colors(k,:),'MarkerSize',10,'parent',axis1); 

    % plot results using outlier removal
    plot(sqrt(Ndofs{2,1}(:,k)),R{2,1}(:,k),'-','DisplayName',dname{1},'linewidth',2,'Color',colors(k,:),'MarkerSize',10,'MarkerFaceColor', 'w','parent',axis2); 
    plot(sqrt(Ndofs{2,2}(:,k)),R{2,2}(:,k),'+','DisplayName',dname{2},'linewidth',2,'Color',colors(k,:),'MarkerSize',10,'MarkerFaceColor', 'w','parent',axis2); 
    plot(sqrt(Ndofs{2,3}(:,k)),R{2,3}(:,k),'--','DisplayName',dname{3},'linewidth',1,'Color',colors(k,:),'MarkerSize',10,'MarkerFaceColor', 'w','parent',axis2);
end
set(axis1,'XScale','log','YScale','log'); box(axis1,'on')
set(axis2,'XScale','log','YScale','log'); box(axis2,'on')

xlim(axis1,[8,64])
xticks(axis1,[8,16,32,64])
ylabel(axis1,'$\| u - u^h \|_{L^2} / \| u \|_{L^2}$ error','FontSize',18,'Interpreter','latex');
ylabel(axis1,'$\| u - u^h \|_{L^2} / \| u \|_{L^2}$','FontSize',18,'Interpreter','latex');
xlabel(axis1,'$\sqrt{N}$','FontSize',18,'Interpreter','latex');
legend1 = legend(axis1,'show');

xlim(axis2,[8,64])
xticks(axis2,[8,16,32,64])
ylabel(axis2,'$\| u - u^h \|_{L^2} / \| u \|_{L^2}$ error','FontSize',18,'Interpreter','latex');
ylabel(axis2,'$\| u - u^h \|_{L^2} / \| u \|_{L^2}$','FontSize',18,'Interpreter','latex');
xlabel(axis2,'$\sqrt{N}$','FontSize',18,'Interpreter','latex');
legend2 = legend(axis2,'show');

%% plot critical time-step

% create figures and axes
fig4  = figure('Position',[0,0,600,400],'Name','Increase in timestep relative to consistent mass (no outlier removal)'); set(gcf, 'Color', 'w');
axis4 = axes('Parent',fig4); hold(axis4,'on');

fig5  = figure('Position',[0,0,600,400],'Name','Increase in timestep relative to consistent mass (using outlier removal)'); set(gcf, 'Color', 'w');
axis5 = axes('Parent',fig5); hold(axis5,'on');


for k=1:length(P)
    p = P(k);
    dname = {sprintf('%s%d', 'Consistent mass p=', p), sprintf('%s%d', 'Approximate mass p=', p), sprintf('%s%d', 'Lumped mass p=', p)};
    
    plot(sqrt(Ndofs{2,2}(:,k)),T{1,2}(:,k)./T{1,1}(:,k),'-','DisplayName',dname{2},'linewidth',2,'Color',colors(k,:),'parent',axis4);
    plot(sqrt(Ndofs{1,1}(:,k)),T{1,3}(:,k)./T{1,1}(:,k),'--','DisplayName',dname{3},'linewidth',1,'Color',colors(k,:),'parent',axis4);

    plot(sqrt(Ndofs{2,2}(:,k)),T{2,2}(:,k)./T{1,2}(:,k),'-','DisplayName',dname{2},'linewidth',2,'Color',colors(k,:),'parent',axis5);
    plot(sqrt(Ndofs{1,1}(:,k)),T{2,3}(:,k)./T{1,3}(:,k),'--','DisplayName',dname{3},'linewidth',1,'Color',colors(k,:),'parent',axis5);
end
set(axis4,'XScale','log','YScale','log'); box(axis4,'on')
set(axis5,'XScale','log','YScale','log'); box(axis5,'on')


xlim(axis4,[8,64]); ylim(axis4, [1,3.0]);
xticks(axis4, [8,16,32,64]); yticks(axis4, 0.0:0.5:5.0);
ylabel(axis4, 'Increase in critical time-step','FontSize',18,'Interpreter','latex');
xlabel(axis4, '$\sqrt{N}$','FontSize',18,'Interpreter','latex');
legend4 = legend(axis4,'show');

xlim(axis5, [8,64]); ylim(axis5, [1,3.0]);
xticks(axis5, [8,16,32,64]); yticks(axis5, 0.0:0.5:5.0);
ylabel(axis5, 'Increase in critical time-step','FontSize',18,'Interpreter','latex');
xlabel(axis5, '$\sqrt{N}$','FontSize',18,'Interpreter','latex');
legend5 = legend(axis5,'show');