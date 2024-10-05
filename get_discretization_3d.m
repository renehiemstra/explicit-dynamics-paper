function [discretization, integrator] = get_discretization_3d(problem, p, ne, Tmax, method, product, flag_mex)
    
    %% Setup discretization

    % pick the discretization method:
    discretization.method = method; % choose: 'standard' or 'dual' or 'lumped'
    discretization.product = product; % choose: 'standard' or 'weighted'

    % construct basis u-dir
    for i = 3:-1:1
        discretization.basis{i}.p = p;
        discretization.basis{i}.nelms = ne(i);
        discretization.basis{i}.kts = knotvector(discretization.basis{i}.p, 0.0 ,1.0 ,ne(i));
        discretization.basis{i}.eval.trialfuns = @(u) trialfuns(discretization.basis{i}, u);
        discretization.basis{i}.eval.testfuns = @(u) testfuns(discretization.basis{i}, u);
        
        % generate quadrature object
        discretization.quadrature.rule{i} = gaussrule(discretization.basis{i});
        m(i) = discretization.quadrature.rule{i}.npoints;
    end
    
    % allocate space for performing sum factorization
    discretization.quadrature.sumfact = {zeros(m), zeros(m), zeros(m), zeros(m)};
    
    % precompute geometric particulars, including first and second derivatives
    discretization.geometry = geometric_particulars_3d(problem.geometry.mapping, discretization.quadrature);

    % pre-compute b-spline basis functions
    for i = 3:-1:1
        discretization.basis{i}.trialfuns = discretization.basis{i}.eval.trialfuns(discretization.quadrature.rule{i}.points);
        discretization.basis{i}.testfuns = discretization.basis{i}.eval.testfuns(discretization.quadrature.rule{i}.points);
        discretization.basis{i}.dim = size(discretization.basis{i}.testfuns{1},2);
        discretization.dims(i) = discretization.basis{i}.dim;
        % 1D mass matrices
        discretization.massmatrix{i} = system_matrix(discretization.basis{i}, discretization.quadrature.rule{i}, 1, 1);
    end 

    % construct mass operators
    % L2 inner product
    if strcmp(discretization.product, 'standard')
        if strcmp(discretization.method, 'standard')
            discretization.mass = assemble_mass_matrix_3d(discretization, flag_mex);
            discretization.chol = chol(discretization.mass);
        elseif strcmp(discretization.method, 'lumped')
            discretization.approxinverse = reshape(1 ./ sum(assemble_mass_matrix_3d(discretization, flag_mex), 2), discretization.basis{1}.dim, discretization.basis{2}.dim, discretization.basis{3}.dim);
        else
            error(strcat('Not implemented for ', discretization.method));
        end
    % weighted L2 inner product
    elseif strcmp(discretization.product, 'weighted')
        if strcmp(discretization.method, 'standard')
            for i = 3:-1:1
                A = discretization.massmatrix{i};
                discretization.basis{i}.scaledtestfuns{1} = discretization.basis{i}.testfuns{1} / A;
                discretization.basis{i}.scaledtestfuns{2} = discretization.basis{i}.testfuns{2} / A;
                discretization.basis{i}.scaledtestfuns{3} = discretization.basis{i}.testfuns{3} / A;
            end
        elseif strcmp(discretization.method, 'lumped')
            for i = 3:-1:1
                A = sum(discretization.massmatrix{i}, 2)';
                discretization.basis{i}.scaledtestfuns{1} = discretization.basis{i}.testfuns{1} ./ A;
                discretization.basis{i}.scaledtestfuns{2} = discretization.basis{i}.testfuns{2} ./ A;
                discretization.basis{i}.scaledtestfuns{3} = discretization.basis{i}.testfuns{3} ./ A;
            end
        elseif strcmp(discretization.method, 'dual')
            for i = 3:-1:1
                A = full(approximate_l2_inverse(discretization.basis{i}.p, discretization.basis{i}.kts, true, true));
                discretization.basis{i}.scaledtestfuns{1} = discretization.basis{i}.testfuns{1} * A;
                discretization.basis{i}.scaledtestfuns{2} = discretization.basis{i}.testfuns{2} * A;
                discretization.basis{i}.scaledtestfuns{3} = discretization.basis{i}.testfuns{3} * A;
            end
        else
            error(strcat('Not implemented for ', discretization.method));
        end
    else
        error('Choose <standard> or <weighted> as discretization.product.');
    end

    %% Setup time integrator 

    % basic Runge-Kutta integrator setup
    integrator = setup_rk_method(discretization);
    
    % final time and safety factor applied to critical time-step
    integrator.tmax = Tmax;
    integrator.safety_factor = 0.5;
    
    % compute the maximum eigenfrequency
    omega_max = estimate_maximum_eigenvalue_3d(discretization, problem, flag_mex);
    
    % compute the critical time-step
    integrator.dtcrit = integrator.cmax / omega_max;
    
    % compute time-step
    nT = ceil(integrator.tmax / (integrator.safety_factor * integrator.dtcrit));
    integrator.dt = integrator.tmax / nT;

end