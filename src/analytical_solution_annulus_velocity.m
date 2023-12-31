function velocity = analytical_solution_annulus_velocity(r,t,theta)
%ANALYTICAL_SOLUTION_ANNULUS_VELOCITY
%    VELOCITY = ANALYTICAL_SOLUTION_ANNULUS_VELOCITY(R,T,THETA)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    15-Apr-2021 01:42:12

velocity = cos(theta.*4.0).*sin(t.*1.106470948850117e+1).*besselj(4,r.*1.106470948850117e+1).*(-1.106470948850117e+1);
