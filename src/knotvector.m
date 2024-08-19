function kts = knotvector(p, a ,b ,ne)
    kts  = [a*ones(1,p) linspace(a, b, ne+1) b*ones(1,p)];
end