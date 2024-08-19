A = rand(2,2);

X = A(:,1);
U = X - [norm(X); 0.0];

I = eye(2);
H = I - (2/dot(U,U)) * (U * U');

R = H * A;
Q = H';

assert(norm(A - Q*R) < 1e-15);