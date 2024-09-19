clear; close; clc;

T = rand(3,4,5);
X_1 = rand(2,3);
X_2 = rand(3,4);
X_3 = rand(4,5);

T = tensorprod(T, X_1, 1, 2);
T = tensorprod(T, X_2, 1, 2);
T = tensorprod(T, X_3, 1, 2);
size(T)