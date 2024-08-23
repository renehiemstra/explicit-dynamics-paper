clear; close; clc;

syms x y t k lambda kappa rho;
assume(x,"real");
assume(y,"real");
assume(t,"real");
assume(kappa,"real");
assume(rho,"real");
assume(k,{'positive','integer'});
assume(lambda,{'positive','integer'});

% displacement, velocity and acceleration
u = sin(k*pi*x) * sin(k*pi*y) * cos(lambda * pi * t);
v = diff(u,t,1);
a = diff(u,t,2);

% compute stiffness
delta_u = divergence(kappa * gradient(u, [x y]), [x y]);

% compute rhs forcing
f = rho * a - delta_u;
simplify(f / u)

