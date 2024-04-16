function [gp] = grevillepts(p, kts)
    n = length(kts) - p - 1;
    gp = zeros(n,1);
    for i=1:n
        gp(i) = sum(kts(i+1:i+p)) / p;
    end
end