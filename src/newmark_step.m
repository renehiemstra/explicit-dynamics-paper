% perform a single Newmark step
function [d_next, v_next, a_next] = newmark_step(d_prev, v_prev, a_prev, t, dt, beta, gamma, solver)
% input:
%   d_prev  -   vector{float64} - displacement previous time step
%   v_prev  -   vector{float64} - velocity previous time step
%   a_prev  -   vector{float64} - acceleration previous time step
%   t       -   float64>0       - current time
%   dt      -   float64>0       - timestep
%   beta    -   0<=float64<=0.5 - Newmark parameter 1
%   gamma   -   0<=float64<=1.0 - Newmark parameter 2
%   solver  -   Function(dF)    - solve for the delta acceleration

    % predictors
    d_next = d_prev + dt*v_prev + (dt^2/2)*(1-2*beta)*a_prev;
    v_next = v_prev + dt*(1-gamma)*a_prev;

    % compute residual forces and delta acceleration
    a_next = solver(d_next, t);

    % update displacement, velocity and acceleration
    v_next = v_next + gamma*dt*a_next;
    d_next = d_next + beta*dt^2*a_next;
    
end



