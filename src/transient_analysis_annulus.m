function [l2e, dt_crit, frame1, frame2, nelms]  = transient_analysis_annulus(D, V, p, nelms, Tmax, tfac, method, remove_outliers, plt_on)
%% discrete solution

outliers = 'false';
if remove_outliers
    outliers = 'true';
end


% radius of 2nd and 4th zero
r_a = 1;
r_b = 17.61596604980483 / 11.064709488501170;

ne_r = nelms;
ne_h = 2*nelms;

% consistent with outliers
F = @(t) 0.0;

% compute extraction operator and mass and stiffness matrices
if strcmp(method, 'standard')
    [kts_r, kts_h, C, L, M, K] = annulus_membrane_model(p, ne_r, ne_h, remove_outliers);
    U = chol(M);
    omega_max = sqrt(abs(eigs(M \ K,1)));
    solver = @(d, v, t) U \ (U' \ (F(t) - K * d));
elseif strcmp(method, 'fast')
    [kts_r, kts_h, C, L, A, K] = annulus_membrane_model_fast(p, ne_r, ne_h, remove_outliers);
    omega_max = sqrt(abs(eigs(A * K,1)));
    solver = @(d, v, t) F(t) - A*(K * d);
elseif strcmp(method, 'lumped')
    [kts_r, kts_h, C, L, A, K] = annulus_membrane_model_lumped(p, ne_r, ne_h, remove_outliers);
    omega_max = sqrt(abs(eigs(A * K,1)));
    solver = @(d, v, t) F(t) - A*(K * d);
end
ndofs_r = length(kts_r)-p-1; 
ndofs_h = length(kts_h)-p-1;
integrate = @(f) test_and_integrate(f, p, kts_r, kts_h, C);
L2_error  = @(u,u_h) l2error_circular_membrane(p, kts_r, kts_h, u, u_h);

% postprocessing
N_r = @(r) spcol(kts_r,p+1,r);
N_h = @(theta) spcol(kts_h,p+1,theta);
Z   = @(d,r,theta) N_r(r) * reshape((C * d(:)), ndofs_h, ndofs_r)' * N_h(theta)';

%% Solve using Newmark method

% final time and timestep as a function of critical time step of
% linear fem with lumped masss
T = Tmax;

% if smaller than 0, take as 90% of critical time step 
if tfac<0
    if p < 3
        time_stepping_method = @(d,v,t,dt,solver) runge_kutta_2_step(d, v, t, dt, solver);
        C_max = 2;
    elseif p<5
        time_stepping_method = @(d,v,t,dt,solver) runge_kutta_4_step(d, v, t, dt, solver);
        C_max = 2.785;
    else
        time_stepping_method = @(d,v,t,dt,solver) runge_kutta_6_step(d, v, t, dt, solver);
        C_max = 3.387;
    end
    dt_crit = C_max / omega_max;
    tfac = 0.5;
end

% dt = tfac * T;
nT = ceil(T / (tfac * dt_crit));
dt    = T / nT;
% beta  = 0.0;
% gamma = 0.5;

% initial condition
d = L \ integrate(@(r,theta) D(r,0,theta)); 
v = L \ integrate(@(r,theta) V(r,0,theta));
nelms = length(d);

% plot initial solution
if plt_on
    % geometry
    X = @(r,theta) r * cos(theta);
    Y = @(r,theta) r * sin(theta);
    R     = linspace(r_a,r_b,50);
    Theta = linspace(0,2*pi,72);
    x = X(R',Theta);
    y = Y(R',Theta);
    z = Z(d,R,Theta);
    
    fig1 = figure(1); set(fig1, 'Color', 'w');
    fig1.Position = [100,100,700,500];
    axis1 = axes('Parent',fig1);
    zmax = max(abs(z(:))) * 1.05;
    contourf(x,y,z,'parent',axis1); 
    hold(axis1,'on');
    plot(x(1,:),y(1,:),'k','linewidth',2,'parent',axis1);
    plot(x(end,:),y(end,:),'k','linewidth',2,'parent',axis1); hold(axis1,'off');
    set(axis1,'Xlim',[-r_b,r_b], 'Ylim', [-1,1]); axis equal; axis off; shading interp; caxis([-zmax zmax]); colormap(jet(20));
    frame1(1) = getframe(fig1);

    fig2 = figure(2); set(fig2, 'Color', 'w');
    fig2.Position = [900,100,700,500];
    axis2 = axes('Parent',fig2);
    contourf(x,y,D(R',0,Theta),'parent',axis2); 
    hold(axis2,'on');
    plot(x(1,:),y(1,:),'k','linewidth',2,'parent',axis2);
    plot(x(end,:),y(end,:),'k','linewidth',2,'parent',axis2); hold(axis2,'off');
    set(axis2,'Xlim',[-r_b,r_b], 'Ylim', [-1,1]); axis equal; axis off; shading interp; caxis([-zmax zmax]); colormap(jet(20));
    frame2(1) = getframe(fig2);
else
    frame1=0;
    frame2=0;
end

%% Newmark scheme

% input:
%   d       -   vector{float64} - initial displacement
%   v       -   vector{float64} - initial velocity
%   a       -   vector{float64} - initial acceleration
%   dt      -   float64>0       - timestep
%   beta    -   0<=float64<=0.5 - Newmark parameter 1
%   gamma   -   0<=float64<=1.0 - Newmark parameter 2
%   solver  -   Function(dF)    - solve for the delta acceleration

i = 2;
for t=dt:dt:T
    
    % status
    fprintf('method = %s, outliers removed = %s, p = %d, n = %d, time-step %d%s%d \n', method{1}, outliers, p, nelms, i, '/', nT);
    
    % compute new displacement, velocity and acceleration
    [d, v] = time_stepping_method(d, v, t, dt, solver);
        
    
    % save figure
    if plt_on
        set(0,'CurrentFigure',fig1);
        z = Z(d,R,Theta);
        contourf(x,y,z,'parent',axis1); 
%         hold(axis1,'on');
%         plot3(x(1,:),y(1,:),z(1,:),'k','linewidth',2,'parent',axis1);
%         plot3(x(end,:),y(end,:),z(1,:),'k','linewidth',2,'parent',axis1); hold(axis1,'off');
        set(axis1,'Xlim',[-r_b,r_b], 'Ylim', [-1,1]); axis equal; axis off; shading interp; caxis([-zmax zmax]); colormap(jet(20));
        frame1(i) = getframe(fig1);
        
        set(0,'CurrentFigure',fig2);
        contourf(x,y,D(R',t,Theta),'parent',axis2); 
        hold(axis2,'on');
        plot(x(1,:),y(1,:),'k','linewidth',2,'parent',axis2);
        plot(x(end,:),y(end,:),'k','linewidth',2,'parent',axis2); hold(axis2,'off');
        set(axis2,'Xlim',[-r_b,r_b], 'Ylim', [-1,1]); axis equal; axis off; shading interp; caxis([-zmax zmax]); colormap(jet(20));
        frame2(i) = getframe(fig2);
    end
    i = i+1;
end

% compute L2 error at time T
l2e = L2_error(@(r,theta) D(r,T,theta), @(r,theta) Z(d,r,theta));