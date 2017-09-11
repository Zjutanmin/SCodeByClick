function [ index_sparsematrix ] = index_eachclass( feature_click, img_label )
%INDEX_EACHCLASS 此处显示有关此函数的摘要
%   此处显示详细说明
    
    [len_query, len_image] = size(feature_click);
%     index_sparsematrix = sparse(len_query,len_query);
    
%     for i = 1:max(len)
%         classindex = find(label == i);
%         num(i) = numel(classindex);
%     end
%     num = arrayfun(@(x) numel(find(label == x)), 1:max(label))';
%     feature_classclick = cell(len_query,1);
    
%     tic
%     for i = 1:len_query
% %         feature_classclick{i} = arrayfun(@(x) sum(feature_click(i,find(label == x))), 1:max(label));
%         for j = 1:max(label)
%             index = find(label == j);
%             feature_classclick{i}(j) = sum(feature_click(i,index));
%         end
%         disp(i);
%     end
%     toc
    
    
    for i = 1:max(img_label)
        index = find(img_label == i);
        feature_classclick(:,i) = sum(feature_click(:,index),2);
        disp(i);
    end
    
    save('feature_classclick.mat','feature_classclick');
    
    f_feature_classclick = full(feature_classclick);
    
%     query_class_count = arrayfun(@(x) numel(find(f_feature_classclick(x,:))),1:len_query);
    each_query_class_index = arrayfun(@(x) find(f_feature_classclick(x,:)),1:len_query,'UniformOutput',false)';
    each_query_class_count = arrayfun(@(x) numel(each_query_class_index{x}),1:len_query)';
    each_query_class_weight = arrayfun(@(x) f_feature_classclick(x,each_query_class_index{x})/sum(f_feature_classclick(x,each_query_class_index{x})),1:len_query,'UniformOutput',false)';
    

    save('each_query_class.mat','each_query_class_count');
    save('each_query_class.mat','each_query_class_index','-append');
    save('each_query_class.mat','each_query_class_weight','-append');
    
    each_class_query_index = arrayfun(@(x) find(f_feature_classclick(:,x)),1:max(img_label),'UniformOutput',false)';
    each_class_query_count = arrayfun(@(x) numel(each_class_query_index{x}),1:max(img_label))';
    
    save('each_class_query.mat','each_class_query_index');
    save('each_class_query.mat','each_class_query_count','-append');
    
    each_class_query_index_e1 = arrayfun(@(x) each_class_query_index{x}(find(each_query_class_count(each_class_query_index{x})==1)),1:max(img_label),'UniformOutput',false)';
    each_class_image_index_e1 = arrayfun(@(x) find(img_label == x),1:max(img_label),'UniformOutput',false)';
    
    each_class_feature_click = cellfun(@(x,y) feature_click(x,y),each_class_query_index_e1,each_class_image_index_e1,'UniformOutput',false);
    
    save('each_class_feature_click.mat','each_class_feature_click');
    save('each_class_feature_click.mat','each_class_query_index_e1','-append');
    save('each_class_feature_click.mat','each_class_image_index_e1','-append');
    
    
    
%     for i = 1:max(label)
%         index_sparsematrix(each_class_query_index{i},each_class_query_index{i}) = 1;
%         index_sparsematrix = triu(index_sparsematrix);
%         disp(i);
%     end


%     [m,n] = cellfun(@test1,each_class_query_index,'UniformOutput',false);

    

%     save('index_sparsematrix.mat','index_sparsematrix');
end


% function [m,n] = test1(x)
%     [m,n] = meshgrid(x,x);
%     [m,n] = deal(reshape(m,[],1),reshape(n,[],1));
% end
