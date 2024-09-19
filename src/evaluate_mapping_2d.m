function [x,y] = evaluate_mapping_2d(map, varargin)

    % compute basis functions
    for i=2:-1:1
        b.u{i} = spcol(map.basis{i}.kts, map.basis{i}.p+1, varargin{i});
    end

    % compute coordinates
    x = tensorconstract(map.cpts{1}, b.u{1}, b.u{2}, 2);
    y = tensorconstract(map.cpts{2}, b.u{1}, b.u{2}, 2);
end