function vol = compute_det(cellmatrix)
    m = size(cellmatrix{1,1});
    vol = zeros(m(1), m(2));
    for l=1:m(2)
        for k=1:m(1)
            A = [cellmatrix{1,1}(k,l) cellmatrix{1,2}(k,l);
                 cellmatrix{2,1}(k,l) cellmatrix{2,2}(k,l)];
            vol(k,l) = det(A);
        end
    end
end