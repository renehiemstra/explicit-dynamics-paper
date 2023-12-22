function [d_next, v_next] = runge_kutta_4_step(d_prev, v_prev, t, dt, solver)
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
d_2 = d_prev + (dt/2)*v_1;
v_2 = v_prev + (dt/2)*a_1;
a_2 = solver(d_2, v_2, t+dt/2);

% third stage
d_3 = d_prev + (dt/2)*v_2;
v_3 = v_prev + (dt/2)*a_2;
a_3 = solver(d_3, v_3, t+dt/2);

% fourth stage
d_4 = d_prev + dt*v_3;
v_4 = v_prev + dt*a_3;
a_4 = solver(d_4, v_4, t+dt);

% compute new displacement and acceleration
d_next = d_prev + (dt/6) * (v_1 + 2*v_2 + 2*v_3 + v_4);
v_next = v_prev + (dt/6) * (a_1 + 2*a_2 + 2*a_3 + a_4);

end