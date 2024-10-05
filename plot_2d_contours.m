function plot_2d_contours(problem, discretization, varargin)

    if nargin==3
        % evaluate the mapping
        gran = varargin{1};
        u = linpts(problem.geometry.mapping.basis{1}, gran(1));
        v = linpts(problem.geometry.mapping.basis{2}, gran(2));
        [x, y] = evaluate_mapping_2d(problem.geometry.mapping, u, v);
        u_bar  = evaluate_field_2d(discretization, u, v);
        % plot contours
        fig = figure; hold on;
        contourf(x, y, u_bar);
        % plot mesh
        % U = unique(discretization.basis{1}.kts);
        % V = unique(discretization.basis{2}.kts);
        % [x, y] = evaluate_mapping_2d(problem.geometry.mapping, U, v);
        % plot(x',y','k','linewidth',2);
        % [x, y] = evaluate_mapping_2d(problem.geometry.mapping, u, V);
        % plot(x,y,'k','linewidth',2);
        axis off;
        export_fig(fig,'solution','-pdf');
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