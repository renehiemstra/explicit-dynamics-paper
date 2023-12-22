function A = approximate_l2_inverse(p, kts, cleft, cright)

    % compute the approximate inverse without boundary enforcement
    A = approxinverse(p+1,p+1,kts(:));

    % specify boundary conditions
    if cleft==true
        A = A - (1/A(1,1)) * A(:,1) * A(1,:);
        A = A(2:end,2:end);
    end
    if cright==true
        A = A - (1/A(end,end)) * A(:,end) * A(end,:);
        A = A(1:end-1,1:end-1);
    end
end


% compute approximate dual of order L for B-spline basis of order m and knot sequence t
function S = approxinverse(L,m,t)

    % initialize
    n = length(t)-m;
    D = speye(n);
    
    % compute S recursively
    S = spdiags(U(m,0,t),0,n,n);
    for v=1:L-1
        D = D * wdifference(t,m+v-1);
        S = S + D * spdiags(U(m,v,t),0,n-v,n-v) * D';
    end
end

% construct n-1 x n sparse representation of scaled difference matrix
function D = wdifference(t, k)
    n = length(t)-k;
    A = k./(t(k+1:end) - t(1:end-k));
    D = spdiags([-A(2:end) A(1:end-1)], -1:0, n, n-1);
end

% compute factors U
function u = U(m,v,t)
    
    % initialize
    n = length(t)-m-v;
    c = factorial(m)*factorial(m-v-1) / (factorial(m+v)*factorial(m+v-1));
    u = zeros(n,1);
    
    % compute factors u depending on multivalued generalized blossoms
    for k=1:n
        u(k) = ((m+v)/(t(k+m+v)-t(k))) * F(v, m+v-1, moment(t(k+1:k+m+v-1), v)) * c;
    end
end

% compute centered moments of degree ?
function sigma = moment(t, v)
    x = mean(t);
    sigma = zeros(2*v,1);
    tau = t-x;
    for l=2:2*v
        sigma(l) = mean(tau.^l);
    end
end

% L² multi-valued generalized blossom labels
function f = F(v, r, sigma)
    if v==0
        f = 1.0;
    elseif v==1
        f = r^2*sigma(2);
    elseif v==2
        f = (1/2)*(r^2*(r^2-3*r+3)*sigma(2)^2 - r^2*(r-1)*sigma(4));
    elseif v==3
        f = (1/6)*(r^3*(r-2)*(r^2-7*r+15)*sigma(2)^3 -3*r^2*(r-2)*(r^2-5*r+10)*sigma(4)*sigma(2) -2*r^2*(3*r^2-15*r+20)*sigma(3)^2 + 2*r^2*(r-1)*(r-2)*sigma(6));
    elseif v==4
        f = (1/24)*(r^4*(r^4 -18*r^3 +125*r^2 -384*r +441)*sigma(2)^4 -6*r^3*(r^4 -16*r^3 +104*r^2 -305*r + 336)*sigma(4)*sigma(2)^2 +3*r^2*(r^4 -14*r^3 +95*r^2 -322*r + 420)*sigma(4)^2 +8*r^2*(r-2)*(r-3)*(r^2-7*r+21)*sigma(6)*sigma(2) -8*r^3*(r-3)*(3*r^2-24*r+56)*sigma(3)^2*sigma(2) +48*r^2*(r-3)*(r^2-7*r+14)*sigma(5)*sigma(3) -6*r^2*(r-1)*(r-2)*(r-3)*sigma(8));
    elseif v==5
        f = (1/120)*(r^5*(r-4)*(r^4-26*r^3+261*r^2-1176*r+2025)*sigma(2)^5 -10*r^4*(r-4)*(r^4-24*r^3+230*r^2-999*r+1674)*sigma(4)*sigma(2)^3 +20*r^3*(r-4)*(r^4-20*r^3+168*r^2-645*r+972)*sigma(6)*sigma(2)^2 +15*r^3*(r-4)*(r^4-22*r^3+211*r^2-942*r+1620)*sigma(4)^2*sigma(2) -20*r^4*(3*r^4-60*r^3+470*r^2-1665*r+2232)*sigma(3)^2*sigma(2)^2 -30*r^2*(r-2)*(r-3)*(r-4)*(r^2-9*r+36)*sigma(8)*sigma(2) -20*r^2*(r-4)*(r^4-18*r^3+173*r^2-828*r+1512)*sigma(6)*sigma(4) +240*r^3*(r^4-19*r^3+143*r^2-493*r+648)*sigma(5)*sigma(3)*sigma(2) +20*r^4*(r-4)*(3*r^2-30*r+83)*sigma(4)*sigma(3)^2 -24*r^2*(5*r^4-90*r^3+655*r^2-2250*r+3024)*sigma(5)^2 -240*r^2*(r-3)*(r-4)*(r^2-9*r+24)*sigma(7)*sigma(3) +24*r^2*(r-1)*(r-2)*(r-3)*(r-4)*sigma(10));
    else
        error('F is not defined for v=$v')
    end
end