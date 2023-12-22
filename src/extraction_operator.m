function C =  extraction_operator(A)
    [n,m] = size(A);
    C = eye(m,m);
    for k=1:n
        Chat = nullspace(A(k,:));
        C = C * Chat;
        A = A * Chat;
    end
end

% ToDo: Add option for two different update algorithms (defining c1)
function C = nullspace(a)

    % obtain permutation for more stable computations
    m = length(a);
    col = zeros(2*m,1);
    row = zeros(2*m,1);
    nzval = zeros(2*m,1);
    [~,p]=sort(abs(a));

    % treat initial zero values
    k=1;j=1;i=1;
    while (j < m+1) && (abs(a(p(i))) < 1e-13)
        row(k) = i; col(k) = j; nzval(k) = 1;
        i=i+1; j=j+1; k=k+1;
    end

    % treat non-zero values
    c2 = 0.0;
    while j < m && (abs(a(p(i+1))) > 1e-13)
        f = a(p(i)) / a(p(i+1));
        %c1 = 1 - c2;
        c1 = 1.0 / (1. + abs(f));
        c2 = -c1 * f;

        row(k) = i; col(k) = j; nzval(k) = c1;
        row(k+1) = i+1; col(k+1) = j; nzval(k+1) = c2;

        i=i+1; j=j+1; k= k+2;
    end

    % treat final zero values
    while i < m && (abs(a(p(i+1))) < 1e-13)
        row(k) = i+1; col(k) = j; nzval(k) = 1;
        i=i+1; j=j+1; k= k+1;
    end

    % create and return sparse matrix
    if (i>m)
        n = m;
    else
        n = m-1;
    end
    C = sparse(p(row(1:k-1)), col(1:k-1), nzval(1:k-1), m, n);
end



















% function C = nullspace(a)
%     m = length(a);
%     C = zeros(m,m-1);
%   
%     i = 1;
%     save = 0.0;
% 
%     for j=1:m-1
%         if a(i)==0
%             C(i,j) = 1;
%         else
%             if a(i+1)==0
%                 i = i+1;
%                 C(i,j) = 1;
%             else
%                 C(i,j)   = 1 - save;
%                 C(i+1,j) = -(a(i)/a(i+1)) *C(i,j);
%                 save = C(i+1,j);
%             end
%         end
%         i = i+1;
%     end
% end