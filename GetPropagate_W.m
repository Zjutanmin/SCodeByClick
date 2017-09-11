function PropagatedVector_W = GetPropagate_W( OriginalVector, alpha, dist, clusterRange, Distype)
    
% % % p_OriginalVector = cellfun(@(x) OriginalVector(:,x(1):x(2)),clusterRange,'UniformOutput',false);
% % % L = cellfun(@(x) x(2)-x(1),clusterRange);
if alpha < 0
    alpha = -alpha;
    N = cellfun(@length, clusterRange, 'UniformOutput', 0);
    clusterratio = cellfun(@(x) alpha*(x-1) /  (alpha*(x-1)+1)*ones(1, x), N, 'UniformOutput', 0);
    temp = bsxfun(@times, OriginalVector, cell2mat(clusterratio));
else
    temp = OriginalVector .* alpha;
end

if nargin < 5
    Distype = 1;
end
esp = 0.001;
switch Distype
    case 1
        dist = exp(-dist);
    case 2
        dist = 1./(dist+esp);
    case 3
        dist = reshape(mapminmax(dist(:),0,1), size(dist));dist = exp(-dist);
    case 4
        dist = reshape(mapminmax(dist(:),0,1), size(dist));dist = 1./(dist+esp);
end

dist = dist - diag(diag(dist));

step = 20;
for j=1:step:size(temp,1)
    k = j:min(j+step-1,size(temp,1));
    for i = 1:length(clusterRange)
        if length(clusterRange{i}) == 1
            continue;
        else
            weight = dist(clusterRange{i},clusterRange{i}) ./ repmat(sum(dist(clusterRange{i},clusterRange{i}),2),1,numel(clusterRange{i}));
            a = repmat(temp(k,clusterRange{i}),1,1,numel(clusterRange{i})) .* permute(repmat(weight,1,1,length(k)),[3,1,2]);
            OriginalVector(k,clusterRange{i}) = OriginalVector(k,clusterRange{i}) + reshape(sum(a,2),size(a,1),size(a,3)) - temp(k,clusterRange{i});
        end
    end
end

PropagatedVector_W = OriginalVector;

end
