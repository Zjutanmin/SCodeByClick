function PropagatedVector_W = GetPropagate_W( OriginalVector, alpha, dist, clusterRange)
    
% % % p_OriginalVector = cellfun(@(x) OriginalVector(:,x(1):x(2)),clusterRange,'UniformOutput',false);
% % % L = cellfun(@(x) x(2)-x(1),clusterRange);

global setting;
temp = OriginalVector .* alpha;
switch setting.type
    case 1
        dist = exp(-dist);
    case 2
        dist = 1./dist;
    case 3
        [m,n] = size(dist);dist = mapminmax(dist(:),0,1);dist = reshape(dist, m,n);dist = exp(-dist);
    case 4
        [m,n] = size(dist);dist = mapminmax(dist(:),0,1);dist = reshape(dist, m,n);dist =  1./dist;
end


dist = dist - diag(ones(size(dist,1),1));
% % % weight = dist ./ repmat(sum(dist,2),1,size(dist,2));

for i = 1:size(temp,1)
% % %     a = weight .* repmat(temp(i,:),size(weight,1),1);
% % %     OriginalVector(i,:) = OriginalVector(i,:) + sum(a,2)';
    for j = 1:length(clusterRange)
        if numel(clusterRange{j}) > 1
            if any(temp(i,clusterRange{j}))
                weight = dist(clusterRange{j},clusterRange{j}) ./ repmat(sum(dist(clusterRange{j},clusterRange{j}),2),1,numel(clusterRange{j}));
                a = weight .* repmat(temp(i,clusterRange{j}),numel(clusterRange{j}),1)';
                OriginalVector(i,clusterRange{j}) = OriginalVector(i,clusterRange{j}) + sum(a,1) - temp(i,clusterRange{j});
            end
        end
    end
end

PropagatedVector_W = OriginalVector;

end