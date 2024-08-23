function plot_2d_contours(problem, discretization, gran)

    % evaluate the mapping
    [x, y] = evaluate_mapping(problem.geometry.mapping, gran);
    u_bar  = evaluate_field(discretization, gran);
   
    % plot mesh
    surf(x, y, u_bar);
end