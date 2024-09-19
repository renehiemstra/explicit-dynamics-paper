function [x,y,z] = evaluate_mapping_3d(map, varargin)

    % compute basis functions
    for i=3:-1:1
        b.u{i} = spcol(map.basis{i}.kts, map.basis{i}.p+1, varargin{i});
    end

    % compute coordinates
    x = tensorconstract(map.cpts{1}, b.u{1}, b.u{2}, b.u{3}, 2);
    y = tensorconstract(map.cpts{2}, b.u{1}, b.u{2}, b.u{3}, 2);
    z = tensorconstract(map.cpts{3}, b.u{1}, b.u{2}, b.u{3}, 2);
end