function X = L2normMatrix(X, tnorm, Normfea)
if nargin < 2
    tnorm = 1;
end
if nargin < 3
    Normfea = 1;
end
if Normfea == 1
idx = find(sum(X.^2, 2));
X(idx,:) = getL2normMatrix(X(idx,:), tnorm);
else if Normfea == -1
        idx = find(sum(X, 2));
        X(idx,:) = getL1normMatrix(X(idx,:), tnorm);
    end
end

function X = getL2normMatrix(X, tnorm)
if issparse(X)
    
    
% [row,col,v] = find(X);
% ind = sub2ind(size(X), row,col);
% Stnorm = sqrt(sum(X.^2, 2));
% X(ind) = v ./ Stnorm(row);


[row,col,v] = find(X);
Stnorm = sqrt(sum(X.^2, 2));
X = sparse(row,col,v ./ Stnorm(row), size(X, 1), size(X, 2));
% nnz(X - Y)
% (ind) = v ./ Stnorm(row);
    
else
X = X ./ repmat(sqrt(sum(X.^2, 2)), [1, size(X, 2)]);

end

if nargin > 1
    X = X ./ sqrt(tnorm);
end

function X = getL1normMatrix(X, tnorm)
if issparse(X)
    
    [row,col,v] = find(X);
% ind = sub2ind(size(X), row,col);
Stnorm = (sum(X, 2));
% X(ind) = v ./ Stnorm(row);



% [row,col,v] = find(X);
% Stnorm = sqrt(sum(X.^2, 2));
X = sparse(row,col,v ./ Stnorm(row), size(X, 1), size(X, 2));


else
X = X ./ repmat((sum(X, 2)), [1, size(X, 2)]);

end

if nargin > 1
    X = X ./ (tnorm);
end