function [X, dX, ddX] = mapping(p, kts, u, xcpts, ycpts)
    % check that input is a rowvector
    assert(size(u,1)==1);

    % evaluation points, including repetitions to compute derivatives of
    % B-splines
    uu = repmat(u,3,1); uu = uu(:);

    % compute B-spline basis functions
    tmp = spcol(kts,p+1,uu); 
    b = tmp(1:3:end,:); db = tmp(2:3:end,:); ddb = tmp(3:3:end,:);

    % plot coordinates
    X = {b * xcpts * b', b * ycpts * b'};

    % compute gradient
    dX = {db * xcpts * b', b * xcpts * db';
          db * ycpts * b', b * ycpts * db'};

    ddX = {ddb * xcpts * b', ddb * ycpts * b';
           b * xcpts * ddb', b * ycpts * ddb';
           db * xcpts * db', db * ycpts * db'};
end