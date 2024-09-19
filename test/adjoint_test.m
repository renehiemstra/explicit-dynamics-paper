clear; close; clc;
A = rand(3,3);

for i=0:2
    for j =0:2
        k1 = mod(i+1, 3);
        k2 = mod(i+2, 3);
        l1 = mod(j+1, 3);
        l2 = mod(j+2, 3);
        Ainv(i+1, j+1) = A(k1+1, l1+1) * A(k2+1, l2+1) - A(k1+1, l2+1) * A(k2+1, l1+1);
    end
end

norm(Ainv - inv(A)' * det(A))
