function T = tensorconstract(T, varargin)
    k = varargin{end};
    for i=1:length(varargin)-1
        T = tensorprod(T, varargin{i}, 1, k);
    end
end