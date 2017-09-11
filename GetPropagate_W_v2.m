function PropagatedVector_W = GetPropagate_W_v2( OriginalVector, alpha, dist, clusterRange, Distype)
    
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
% % % weight = dist ./ repmat(sum(dist,2),1,size(dist,2));
% if isempty(gcp('nocreate'))
%     parpool('local',4);
% end
t = [];
% par
for i = 1:size(temp,1)
% % %     a = weight .* repmat(temp(i,:),size(weight,1),1);
% % %     OriginalVector(i,:) = OriginalVector(i,:) + sum(a,2)';
%     for j = 1:length(clusterRange)
%         if numel(clusterRange{j}) > 1
%             if any(temp(i,clusterRange{j}))
%                 weight = dist(clusterRange{j},clusterRange{j}) ./ repmat(sum(dist(clusterRange{j},clusterRange{j}),2),1,numel(clusterRange{j}));
%                 a = weight .* repmat(temp(i,clusterRange{j}),numel(clusterRange{j}),1)';
%                 OriginalVector(i,clusterRange{j}) = OriginalVector(i,clusterRange{j}) + sum(a,1) - temp(i,clusterRange{j});
%             end
%         end
%     end
    disp(i);
%i = 1419;   
    t(i,:) = propagate(clusterRange,dist,temp(i,:), OriginalVector(i,:));
end

% delete(gcp);

% PropagatedVector_W = OriginalVector;
PropagatedVector_W = t;

end


function vector = propagate(clusterRange,dist,temp, OriginalVector)
% if nnz(temp)==0
%     vector = OriginalVector;
% else
    for j = 1:length(clusterRange)
        if numel(clusterRange{j}) > 1
            if any(temp(clusterRange{j}))
                weight = dist(clusterRange{j},clusterRange{j}) ./ repmat(sum(dist(clusterRange{j},clusterRange{j}),2),1,numel(clusterRange{j}));
                a = weight .* repmat(temp(clusterRange{j}),numel(clusterRange{j}),1)';
                vector(clusterRange{j}) = OriginalVector(clusterRange{j}) + sum(a,1) - temp(clusterRange{j});
            else
                vector(clusterRange{j}) = 0;
            end
        end
    end
    
end
% end