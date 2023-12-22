function [d_next, v_next] = runge_kutta_6_step(d_prev, v_prev, t, dt, solver)
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
d_2 = d_prev + (dt/4)*v_1;
v_2 = v_prev + (dt/4)*a_1;
a_2 = solver(d_2, v_2, t+dt/4);

% third stage
d_3 = d_prev + (dt/8)*(v_1+v_2);
v_3 = v_prev + (dt/8)*(a_1+a_2);
a_3 = solver(d_3, v_3, t+dt/4);

% fourth stage
d_4 = d_prev + (dt/2)*(-v_2 + 2*v_3);
v_4 = v_prev + (dt/2)*(-a_2 + 2*a_3);
a_4 = solver(d_4, v_4, t+dt/2);

% fifth stage
d_5 = d_prev + (dt/16) * (3*v_1 + 9*v_4);
v_5 = v_prev + (dt/16) * (3*a_1 + 9*a_4);
a_5 = solver(d_5, v_5, t+3*dt/4);

% sixth stage
d_6 = d_prev + (dt/7) * (-3*v_1 + 2*v_2 + 12*v_3 - 12*v_4 + 8*v_5);
v_6 = v_prev + (dt/7) * (-3*a_1 + 2*a_2 + 12*a_3 - 12*a_4 + 8*a_5);
a_6 = solver(d_6, v_6, t+dt);

% compute new displacement and acceleration
d_next = d_prev + (dt/90) * (7*v_1 + 32*v_3 + 12*v_4 + 32*v_5 + 7*v_6);
v_next = v_prev + (dt/90) * (7*a_1 + 32*a_3 + 12*a_4 + 32*a_5 + 7*a_6);

end