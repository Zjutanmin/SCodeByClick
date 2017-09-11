function [cluster_idx] = ImageClustering(Keachone_img)
      
%     for i = 1:max(img_label) 
%         index = find(data_label == i);
%         fea = data_fea(index,:)';
%         [C, A] = vl_kmeans(fea,Keachone(i));
%     end

    load('image_click_Dog283_0_CNN_Alex1_ND_S_S1_data_normal.mat');
% % %         data_fea = data_fea ./ repmat(sqrt(sum(data_fea.^2,2)),1,size(data_fea,2));
% % %         dist = EuDist2(data_fea,data_fea);
    
    [~, A] = arrayfun(@(x) vl_kmeans(data_fea(data_label==x,:)',Keachone_img(x)),1:max(data_label),'UniformOutput',false);
%     [~, A] = arrayfun(@(x) kmeans(data_fea(find(data_label==x),:)',Keachone_img(x)),1:max(data_label),'UniformOutput',false);
    cluster_idx = arrayfun(@(x) A{x}+sum(Keachone_img(1:x-1)),1:max(data_label),'UniformOutput',false)';
    cluster_idx = cell2mat(cluster_idx')';
end