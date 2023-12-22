function [d_next, v_next] = runge_kutta_2_step(d_prev, v_prev, t, dt, solver)
% input:
%   d_prev  -   vector{float64} - displacement previous time step
%   v_prev  -   vector{float64} - velocity previous time step
%   t       -   float64>0       - current time
%   dt      -   float64>0       - timestep
%   solver  -   Function(dF)    - solve for the delta acceleration

% first stage
d_1 = d_prev;
v_1 = v_prev;
a_1 = solver(d_1, v_1, t);

% second stage
d_2 = d_1 + (dt/2)*v_1;
v_2 = v_1 + (dt/2)*a_1;
a_2 = solver(d_2, v_2, t+dt/2);

% compute new displacement and acceleration
d_next = d_prev + dt*v_2;
v_next = v_prev + dt*a_2;

end