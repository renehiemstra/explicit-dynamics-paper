function [x,y] = evaluate_mapping(map, granularity)

    % initialize
    u.p = map.basis{1}.p; u.kts = map.basis{1}.kts;
    v.p = map.basis{2}.p; v.kts = map.basis{2}.kts;

    % compute basis functions
    b.u = spcol(u.kts, u.p+1, linspace(u.kts(1), u.kts(end), granularity(1)));
    b.v = spcol(v.kts, v.p+1, linspace(v.kts(1), v.kts(end), granularity(2)));

    % compute coordinates
    x = b.u * map.cpts{1} * b.v';
    y = b.u * map.cpts{2} * b.v';
end