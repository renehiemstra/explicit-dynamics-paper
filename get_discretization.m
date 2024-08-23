function [discretization, integrator] = get_discretization(problem, p, ne, Tmax, method, product, flag_mex)
    
    %% Setup discretization

    % pick the discretization method:
    discretization.method = method; % choose: 'standard' or 'dual' or 'lumped'
    discretization.product = product; % choose: 'standard' or 'weighted'

    % construct basis u-dir
    discretization.basis{1}.p = p;
    discretization.basis{1}.nelms = ne(1);
    discretization.basis{1}.kts = knotvector(discretization.basis{1}.p, 0.0 ,1.0 ,ne(1));
    discretization.basis{1}.eval.trialfuns = @(u) trialfuns(discretization.basis{1}, u);
    discretization.basis{1}.eval.testfuns = @(u) testfuns(discretization.basis{1}, u);
    
    % construct basis v-dir
    discretization.basis{2}.p = p;
    discretization.basis{2}.nelms = ne(2);
    discretization.basis{2}.kts = knotvector(discretization.basis{2}.p, 0.0 ,1.0 ,ne(2));
    discretization.basis{2}.eval.trialfuns = @(v) trialfuns(discretization.basis{2}, v);
    discretization.basis{2}.eval.testfuns = @(v) testfuns(discretization.basis{2}, v);
    
    % generate quadrature object
    discretization.quadrature.rule{1} = gaussrule(discretization.basis{1});
    discretization.quadrature.rule{2} = gaussrule(discretization.basis{2});
    m = [discretization.quadrature.rule{1}.npoints, discretization.quadrature.rule{2}.npoints];
    % allocate space for performing sum factorization
    discretization.quadrature.sumfact = {zeros(m(1), m(2)), zeros(m(1), m(2)), zeros(m(1), m(2))};
    
    % precompute geometric particulars, including first and second derivatives
    discretization.geometry = geometric_particulars(problem.geometry.mapping, discretization.quadrature);
    
    % pre-compute b-spline basis functions
    % basis in u-dir
    discretization.basis{1}.trialfuns = discretization.basis{1}.eval.trialfuns(discretization.quadrature.rule{1}.points);
    discretization.basis{1}.testfuns = discretization.basis{1}.eval.testfuns(discretization.quadrature.rule{1}.points);
    discretization.basis{1}.dim = size(discretization.basis{1}.testfuns{1},2);
    % basis in v-dir
    discretization.basis{2}.trialfuns = discretization.basis{2}.eval.trialfuns(discretization.quadrature.rule{2}.points);
    discretization.basis{2}.testfuns = discretization.basis{2}.eval.testfuns(discretization.quadrature.rule{2}.points);
    discretization.basis{2}.dim = size(discretization.basis{2}.testfuns{1},2);
    
    % compute the mass matrix
    discretization.massmatrix{1} = system_matrix(discretization.basis{1}, discretization.quadrature.rule{1}, 1, 1);
    discretization.massmatrix{2} = system_matrix(discretization.basis{2}, discretization.quadrature.rule{2}, 1, 1);
    
    if strcmp(discretization.product, 'standard')
        if strcmp(discretization.method, 'standard')
            discretization.mass = assemble_mass_matrix(discretization, flag_mex);
            discretization.chol = chol(discretization.mass);
        elseif strcmp(discretization.method, 'lumped')
            discretization.approxinverse = reshape(1 ./ sum(assemble_mass_matrix(discretization, flag_mex), 2), discretization.basis{1}.dim, discretization.basis{2}.dim);
        else
            error(strcat('Not implemented for ', discretization.method));
        end
    elseif strcmp(discretization.product, 'weighted')
        if strcmp(discretization.method, 'standard')
            discretization.basis{1}.testfuns = discretization.basis{1}.testfuns / discretization.massmatrix{1};
            discretization.basis{2}.testfuns = discretization.basis{2}.testfuns / discretization.massmatrix{2};
        elseif strcmp(discretization.method, 'lumped')
            discretization.basis{1}.testfuns = discretization.basis{1}.testfuns ./ sum(discretization.massmatrix{1}, 2)';
            discretization.basis{2}.testfuns = discretization.basis{2}.testfuns ./ sum(discretization.massmatrix{2}, 2)';
        elseif strcmp(discretization.method, 'dual')
            A_1 = full(approximate_l2_inverse(discretization.basis{1}.p, discretization.basis{1}.kts, true, true));
            A_2 = full(approximate_l2_inverse(discretization.basis{2}.p, discretization.basis{2}.kts, true, true));
            discretization.basis{1}.scaledtestfuns{1} = discretization.basis{1}.testfuns{1} * A_1;
            discretization.basis{1}.scaledtestfuns{2} = discretization.basis{1}.testfuns{2} * A_1;
            discretization.basis{1}.scaledtestfuns{3} = discretization.basis{1}.testfuns{3} * A_1;     
            discretization.basis{2}.scaledtestfuns{1} = discretization.basis{2}.testfuns{1} * A_2;
            discretization.basis{2}.scaledtestfuns{2} = discretization.basis{2}.testfuns{2} * A_2;
            discretization.basis{2}.scaledtestfuns{3} = discretization.basis{2}.testfuns{3} * A_2;
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
    omega_max = estimate_maximum_eigenvalue(discretization, problem);
    
    % compute the critical time-step
    integrator.dtcrit = integrator.cmax / omega_max;
    
    % compute time-step
    nT = ceil(integrator.tmax / (integrator.safety_factor * integrator.dtcrit));
    integrator.dt = integrator.tmax / nT;

end