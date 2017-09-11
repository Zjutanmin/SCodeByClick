function  [QueryFea] = getqueryfeature(sparseOrg, cluster_idx, type,alpha,dist)
global setting;
% cluster_idx = cluster_idx(1:390);
    switch type
        case 'O'   
            QueryFea = sparseOrg;
            issparse(QueryFea)
            
        case 'C'   
            QueryFea = arrayfun(@(x) sum(sparseOrg(:,cluster_idx==x),2),min(cluster_idx):max(cluster_idx),'UniformOutput',false);
            QueryFea = cell2mat(QueryFea);
            issparse(QueryFea)
        case 'P'   % %��ֵ����
            dist = ones([size(sparseOrg, 2), size(sparseOrg, 2)]);
            clusterRange = arrayfun(@(x) [find(cluster_idx==x)],min(cluster_idx):max(cluster_idx),'UniformOutput',false);
            QueryFea = GetPropagate_W(full(sparseOrg),alpha,dist,clusterRange, setting.type);
% %             clusterRange = arrayfun(@(x) [find(cluster_idx==x,1,'first'),find(cluster_idx==x,1,'last')],min(cluster_idx):max(cluster_idx),'UniformOutput',false);
% %             t = cellfun(@(x) numel(x)==0,clusterRange);
% %             t = find(t == 1);
% %             for i = 1:length(t)
% %                 clusterRange{t(i)} = [1,0];
% %             end
% %             Original_index = arrayfun(@(x) find(cluster_idx == x),min(cluster_idx):max(cluster_idx),'UniformOutput',false);
% %             Original_index = cell2mat(Original_index')';
% %             QueryFea = GetPropagateClick(sparseOrg,alpha,clusterRange,Original_index);
% %             issparse(QueryFea)
            
        case 'WP'     % %Ȩ�ش���
%             clusterRange = arrayfun(@(x) [find(cluster_idx==x,1,'first'),find(cluster_idx==x,1,'last')],min(cluster_idx):max(cluster_idx),'UniformOutput',false);
%             QueryFea = GetPropagate_W(sparseOrg,alpha,dist,clusterRange);
            
            clusterRange = arrayfun(@(x) [find(cluster_idx==x)],min(cluster_idx):max(cluster_idx),'UniformOutput',false);
            QueryFea = GetPropagate_W(full(sparseOrg),alpha,dist,clusterRange, setting.type);
    end
    
    if length(find(isnan(QueryFea))) > 1
        error('Data error\n')
    end
    
    

end

