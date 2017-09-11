% % % 计算同一类英文和非英文样本的距离平均值


load('nonE_query_idx.mat');
load('result/N-C-k-10-20-cluster1.mat');
nonE_query_cluster = unique(IDX(nonE_query_idx));
load('ImgC/ImgC-N-10.mat');
load('img_Fea_Clickcount');
% % % image_idx = arrayfun(@(x) find(cluster_idx==x),min(cluster_idx):max(cluster_idx),'UniformOutput',false);


% % % get the matrix of click vector
clickvector = arrayfun(@(x) sum(img_Fea_Clickcount(cluster_idx==x,:),1)',min(cluster_idx):max(cluster_idx),'UniformOutput',false);
clickvector = cell2mat(clickvector);
clear img_Fea_Clickcount;


meandist = zeros(length(nonE_query_cluster),1);
meandist_normal = zeros(length(nonE_query_cluster),1);
A = cell(length(nonE_query_cluster),1);
B = cell(length(nonE_query_cluster),1);
A_normal = cell(length(nonE_query_cluster),1);
B_normal = cell(length(nonE_query_cluster),1);
for i = 1:length(nonE_query_cluster)
    q = find(IDX==nonE_query_cluster(i));
    nonE_query = nonE_query_idx(find(IDX(nonE_query_idx)==nonE_query_cluster(i)));
    E_query = setdiff(q,nonE_query);
    
    A{i,1} = clickvector(nonE_query,:);
    B{i,1} = clickvector(E_query,:);
    
    meandist(i,1) = mean(mean(EuDist2(A{i,1},B{i,1})));
    
    A_normal{i,1} = A{i,1}./repmat(sqrt(sum(A{i,1}.^2,2)),1,size(A{i,1},2));
    B_normal{i,1} = B{i,1}./repmat(sqrt(sum(B{i,1}.^2,2)),1,size(B{i,1},2));
    
    meandist_normal(i,1) = mean(mean(EuDist2(A_normal{i,1},B_normal{i,1})));
    
    disp(i);
end
save('meandist.mat','meandist');
save('meandist_normal','meandist_normal');
save('clickvector','A','B');
save('clickvector_normal','A_normal','B_normal');
clear IDX;


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 



load('nonE_query_idx.mat');
load('result/N-C-k-10-20-cluster1-0.mat');
nonE_query_cluster = unique(IDX(nonE_query_idx));

clickvector = [];
sum = 0;
for i = 1:283
    load(['1-283/' num2str(i) '.mat']);
    sum = sum + length(clus);
    disp(i);
end

meandist_text = zeros(length(nonE_query_cluster),1);
meandist_normal_text = zeros(length(nonE_query_cluster),1);
A = cell(length(nonE_query_cluster),1);
B = cell(length(nonE_query_cluster),1);
A_normal = cell(length(nonE_query_cluster),1);
B_normal = cell(length(nonE_query_cluster),1);
for i = 1:length(nonE_query_cluster)
    q = find(IDX==nonE_query_cluster(i));
    nonE_query = nonE_query_idx(find(IDX(nonE_query_idx)==nonE_query_cluster(i)));
    E_query = setdiff(q,nonE_query);
    
    A{i,1} = clickvector(nonE_query,:);
    B{i,1} = clickvector(E_query,:);
    
    meandist_text(i,1) = mean(mean(EuDist2(A{i,1},B{i,1})));
    
    A_normal{i,1} = A{i,1}./repmat(sqrt(sum(A{i,1}.^2,2)),1,size(A{i,1},2));
    B_normal{i,1} = B{i,1}./repmat(sqrt(sum(B{i,1}.^2,2)),1,size(B{i,1},2));
    
    meandist_normal_text(i,1) = mean(mean(EuDist2(A_normal{i,1},B_normal{i,1})));
    
    disp(i);
end
