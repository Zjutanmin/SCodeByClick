%求每个query点击数最大的一个图像小类的索引


method = 'N';
k_cluster_I = 20;

load('each_class_feature_click.mat');
Neachone_img = cellfun(@(x) numel(x),each_class_image_index_e1);

Keachone_img = getKeachclass_img(method, k_cluster_I, Neachone_img);
% % % load('image_click_Dog283_0_CNN_Alex1_ND_S_S1_data.mat');
% % % [~, A] = arrayfun(@(x) vl_kmeans(data_fea(data_label==x,:)',Keachone_img(x)),1:max(data_label),'UniformOutput',false);
% % % A = arrayfun(@(x) A{x}+sum(Keachone_img(1:x-1)),1:max(data_label),'UniformOutput',false)';
% % % % % % for i = 1:length(each_class_query_index_e1)
% % %     


sname = ['ImgC/ImgC-' method '-' num2str(k_cluster_I) '.mat'];
load(sname);
load('img_Fea_Clickcount.mat');

feature_click = img_Fea_Clickcount';

for i = 1:max(cluster_idx)
    index = find(cluster_idx == i);
    feature_classclick(:,i) = sum(feature_click(:,index),2);
    disp(i);
end

[~, cluster_max_idx] = max(feature_classclick,[],2);

cluster_max_idx = cellfun(@(x) cluster_max_idx(x),each_class_query_index_e1,'UniformOutput',false);
save('cluster_max_idx.mat','cluster_max_idx');


% % % % for i = 1:length(each_class_query_index_e1)
% % % %     each_query_lessclass_index  = arrayfun(@(x) find(feature_classclick(x,:)),each_class_query_index_e1{i},'UniformOutput',false)';
% % % %     each_query_lessclass_count = arrayfun(@(x) numel(each_query_lessclass_index{x}),1:length(each_class_query_index_e1{i}))';
% % % %     arrayfun(@(x) each_query_lessclass_index(each_query_lessclass_count>1))
% % % %     
% % % %     each_lessclass_query_index = arrayfun(@(x) find(feature_classclick(:,x)),1:max(cluster_idx),'UniformOutput',false)';
% % % %     each_lessclass_query_count = arrayfun(@(x) numel(each_lessclass_query_index{x}),1:max(cluster_idx))';
% % % %     
% % % % end
% % % % 
% % % % 
% % % % % eachless_query_class_index = arrayfun(@(x) find(f_feature_classclick(x,:)),1:len_query,'UniformOutput',false)';
% % % % % eachless_query_class_count = arrayfun(@(x) numel(eachless_query_class_index{x}),1:len_query)';
% % % % % eachless_query_class_weight = arrayfun(@(x) f_feature_classclick(x,eachless_query_class_index{x})/sum(f_feature_classclick(x,eachless_query_class_index{x})),1:len_query,'UniformOutput',false)';
% % % % %     