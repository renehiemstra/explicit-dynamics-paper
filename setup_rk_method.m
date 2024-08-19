function rk = setup_rk_method(discretization)
    p = max(discretization.basis{1}.p, discretization.basis{2}.p);
    if (p < 3)
        rk.solver = @(d,v,t,dt,solver) runge_kutta_2_step(d, v, t, dt, solver);
        rk.cmax = 2;
    elseif p<5
        rk.solver = @(d,v,t,dt,solver) runge_kutta_4_step(d, v, t, dt, solver);
        rk.cmax = 2.785;
    else
        rk.solver = @(d,v,t,dt,solver) runge_kutta_6_step(d, v, t, dt, solver);
        rk.cmax = 3.387;
    end
end