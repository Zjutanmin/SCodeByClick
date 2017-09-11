function x =getSim(x)
    if ~isempty(x)
        x = exp(-EuDist2(x,x));
        eps = 1e-3;
        x = reshape(mapminmax(x(:)', eps, 1),size(x));
    end
end