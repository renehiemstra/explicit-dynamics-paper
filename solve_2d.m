function discretization = solve_2d(problem, discretization, integrator, flag_mex)

    for i=2:-1:1
        P(i) = discretization.basis{i}.p;
    end
    p = max(P);

    % get initial quess using L2-projection
    d = project_2d(discretization, @(x,y) problem.solution.displacement(x,y,0.0));
    v = project_2d(discretization, @(x,y) problem.solution.velocity(x,y,0.0));
    
    % choose solver object
    if strcmp(discretization.product, 'standard')
        % standard - L2 inner product
        solver = @(d, v, t) assemble_standard_forcing_vector_2d(d, discretization, problem.material, @(x, y) problem.solution.forcing(x, y, t), flag_mex);
    elseif strcmp(discretization.product, 'weighted')
        % weighted L2 inner product
        solver = @(d, v, t) assemble_weighted_forcing_vector_2d(d, discretization, problem.material, @(x, y) problem.solution.forcing(x, y, t), flag_mex);
    end
    

    % Runge-Kutta time integration
    T = integrator.tmax;
    dt = integrator.dt;
    i = 1;
    for t=0:dt:T-dt
        
        % status
        if mod(i, 100)==0
            fprintf('time integrator, p = %d, mesh = %d x %d, time-step %d%s%d \n', p, discretization.basis{1}.nelms, discretization.basis{2}.nelms,  i, '/', integrator.tmax/integrator.dt);
        end

        % compute new displacement, velocity and acceleration
        [d, v] = integrator.solver(d, v, t, dt, solver);
    
        % increase counter
        i = i+1;
    end

    discretization.coeffs = d;
end