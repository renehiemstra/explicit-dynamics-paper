function plot_2d_contours(problem, discretization, varargin)

    if nargin==3
        % evaluate the mapping
        gran = varargin{1};
        u = linpts(problem.geometry.mapping.basis{1}, gran(1));
        v = linpts(problem.geometry.mapping.basis{2}, gran(2));
        [x, y] = evaluate_mapping_2d(problem.geometry.mapping, u, v);
        u_bar  = evaluate_field_2d(discretization, u, v);
       
        % plot mesh
        surf(x, y, u_bar);
    elseif nargin==4
        gran = varargin{1};
        crossection = varargin{2};
        u = linpts(problem.geometry.mapping.basis{1}, gran(1));
        v = linpts(problem.geometry.mapping.basis{2}, gran(2));
        w = [crossection];
        % evaluate the mapping
        [x, y, z] = evaluate_mapping_3d(problem.geometry.mapping, u, v, w);
        u_bar  = evaluate_field_3d(discretization, u, v, w);
        % plot mesh
        surf(x, y, u_bar);
    end
end

function u = linpts(basis, gran)
    u = linspace(basis.kts(1), basis.kts(end), gran);
end