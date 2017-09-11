function IDX = mergeQ(sname, method,theta,alpha,P,Distype,cluster_method, ...
    type, T,gamma, textweight,normweight,k_cluster_Q)               
if nargin < 11
    textweight = 1;
end
if nargin < 12
    normweight = 1;
end
if nargin < 13
    k_cluster_Q = 20;
end
k_cluster_I = 20;

type_empty_cluster = 1;
global setting;
setting.type = Distype;


try
    load(fullfile('result', [sname, 'c-', num2str(type_empty_cluster), '-IDX.mat']), 'IDX');
catch
    
    addpath(genpath('ksvd+omp'))
    addpath('data')
    addpath('TextFea/image_click_Dog283_0_QueryFea')
    addpath(genpath('vlfeat'))
    vl_setup
    
    h = tic;
    load('img_Fea_Clickcount.mat');
    load('each_class_feature_click.mat');
    load('each_class_query.mat');
    load('each_query_class.mat');
    feature_click = img_Fea_Clickcount';
    Neachone_img = cellfun(@(x) numel(x),each_class_image_index_e1);
    Neachone_query = cellfun(@(x) numel(x),each_class_query_index_e1);
        
    
    each_class_index_query_g1_empty = find(Neachone_query == 0);
    len_empty = numel(each_class_index_query_g1_empty);
    each_class_query_g1_empty_idx = arrayfun(@(x) each_class_query_index{x},each_class_index_query_g1_empty,'UniformOutput',false);
    query_g1_empty_idx = cell2mat(each_class_query_g1_empty_idx);
    query_g1_empty_idx = unique(query_g1_empty_idx);
    
%     query_g1_nonem_idx = find(each_query_class_count>1);
%     query_g1_nonem_idx = query_g1_nonem_idx(find(ismember(query_g1_nonem_idx,query_g1_empty_idx)==0));

    
    h1= toc(h);
    sprintf('load %d', h1)
    whos
    h = tic;
    
    try
        load(['ImgC/' 'ImgC-', method, '-', num2str(k_cluster_I), '.mat'], 'cluster_idx')
    catch
        Keachone_img = getKeachclass_img(method, k_cluster_I, Neachone_img);
        cluster_idx = ImageClustering(Keachone_img);
        save(['ImgC/' 'ImgC-', method, '-', num2str(k_cluster_I), '.mat'], 'cluster_idx')
        disp('save successful');
    end
    
    h1= toc(h);
    sprintf('ImagClustering %d', h1)
    whos
    h = tic;
    

    [sort_cluster_idx, sort_idx] = sort(cluster_idx);
    sort_feature_click = feature_click(:,sort_idx);
    
    
                
%                 h1= toc(h);
%     sprintf('getqueryfeature %d', h1)
%     whos
%     h = tic;
    
        %Neachone_query = query number in each category
    Keachone_query = getKeachclass_query(method, k_cluster_Q, Neachone_query);
    
    if type == 'WP'
        load('image_click_Dog283_0_CNN_Alex1_ND_S_S1_data_normal_label.mat');
        load('image_click_Dog283_0_CNN_Alex1_ND_S_S1_data_normal_fea1.mat');
        load('image_click_Dog283_0_CNN_Alex1_ND_S_S1_data_normal_fea2.mat');
        load('image_click_Dog283_0_CNN_Alex1_ND_S_S1_data_normal_fea3.mat');
        load('image_click_Dog283_0_CNN_Alex1_ND_S_S1_data_normal_fea4.mat');
        load('image_click_Dog283_0_CNN_Alex1_ND_S_S1_data_normal_fea5.mat');
        load('image_click_Dog283_0_CNN_Alex1_ND_S_S1_data_normal_fea6.mat');
        load('image_click_Dog283_0_CNN_Alex1_ND_S_S1_data_normal_fea7.mat');
        load('image_click_Dog283_0_CNN_Alex1_ND_S_S1_data_normal_fea8.mat');
        load('image_click_Dog283_0_CNN_Alex1_ND_S_S1_data_normal_fea9.mat');
        load('image_click_Dog283_0_CNN_Alex1_ND_S_S1_data_normal_fea10.mat');
        load('image_click_Dog283_0_CNN_Alex1_ND_S_S1_data_normal_fea11.mat');
        data_fea = [ data_fea1,data_fea2,data_fea3,data_fea4,data_fea5,data_fea6,data_fea7,data_fea8,data_fea9,data_fea10,data_fea11];
        clear data_fea1 data_fea2 data_fea3 data_fea4 data_fea5 data_fea6 data_fea7 data_fea8 data_fea9 data_fea10 data_fea11
        data_fea = full(data_fea);
% % %         data_fea = data_fea ./ repmat(sqrt(sum(data_fea.^2,2)),1,size(data_fea,2));
% % %         dist = EuDist2(data_fea,data_fea);
    end
    
   

    try 
        load(fullfile('result', [sname, '-KNN.mat']),'KNNIDX','KNNDis');
    catch
        KNNIDX = ones(length(each_query_class_index), 1) * 0;
        KNNDis = ones(length(each_query_class_index), 1) * inf;
        for i = 1:length( Keachone_query)
            h = tic;
            if ~isempty(each_class_query_index_e1{i})
                Q_index_E1 = each_class_query_index_e1{i};
                Q_index_E2 = setdiff(setdiff(each_class_query_index{i}', query_g1_empty_idx'), Q_index_E1);
                Q_index_E = [Q_index_E1; Q_index_E2'];
                range1 = 1:length(Q_index_E1);
                range2 = length(Q_index_E1)+1:length(Q_index_E);
                if type == 'WP'
                    t_dist = EuDist2(data_fea(data_label==i,:),data_fea(data_label==i,:));
                else
                    t_dist = [];
                end
                %%%pro
                each_class_idx = sum(Neachone_img(1:i))-Neachone_img(i)+1:sum(Neachone_img(1:i));
% %                 %%%normal
% %                 QueryFea = sort_feature_click(Q_index_E,each_class_idx) ./ repmat(sqrt(sum(sort_feature_click(Q_index_E,each_class_idx).^2,2)),1,size(sort_feature_click(Q_index_E,each_class_idx),2));
%                 tic
%                 QueryFea = getqueryfeature(QueryFea, sort_cluster_idx(each_class_idx), type,alpha,t_dist,normweight);
%                 toc
                ClickCount = sum(sort_feature_click(Q_index_E,each_class_idx),2);
                QueryFea = getqueryfeature(sort_feature_click(Q_index_E,each_class_idx), sort_cluster_idx(each_class_idx), ...
                    type,alpha,t_dist);
                %%clustering
                temp_QueryFea = QueryFea(range1,:)';ClickCount = ClickCount(range1);
                [cluster_vector,cluster_idx_q] = Myclustering(temp_QueryFea,...
                    Keachone_query,i,cluster_method,gamma,P,theta,T, textweight,normweight,ClickCount);
                
                %%[�ֵ䣬����coding]
                %         [cluster_idx_q, cluster_vector] = kmeans(QueryFea(range1,:),Keachone_query(i));
                cluster_idx_q = cluster_idx_q + max(KNNIDX);
                
                KNNIDX(Q_index_E(range1)) = cluster_idx_q;
                %             cluster_vector = cluster_vector(len_clickvector,:);


                %%1nn
                range = Q_index_E(range2);
                [dist, idx] = min(EuDist2(QueryFea(range2,:), cluster_vector'), [], 2);
                idx = double(cluster_idx_q(idx'))';
                tmp= KNNDis(range) < dist;
                KNNIDX(range) =  (tmp) .* KNNIDX(range) + ~tmp .* idx;
                %         KNNDis(range2) =  (tmp) .* KNNDis(range2) + ~tmp .* dist;
                KNNDis(range(tmp==0)) = dist(tmp==0);
                %             IDX(Q_index_E) = KNNIDX;
                disp(i);
            end
            h1 = toc(h);
            sprintf('%d : %d',i,h1);
        end
        h1= toc(h);
        sprintf(' QueryClustering %d', h1)
        whos
        h = tic;
    
        save(fullfile('result', [sname, '-KNN.mat']),'KNNIDX','KNNDis');
        disp('save successful');
    end
    
    
    if type_empty_cluster == 1
        class_idx = arrayfun(@(x) sum(Neachone_img(1:x))-Neachone_img(x)+1:sum(Neachone_img(1:x)),each_class_index_query_g1_empty,'UniformOutput',false);
        class_idx = cell2mat(class_idx')';
        
        if type == 'WP'
            t_dist = EuDist2(data_fea(class_idx,:),data_fea(class_idx,:));
        else
            t_dist = [];
        end
        
%         %%%normal
%         temp_QueryFea = sort_feature_click(query_g1_empty_idx,class_idx)./repmat(sqrt(sum(sort_feature_click(query_g1_empty_idx,class_idx).^2,2)),1,size(sort_feature_click(query_g1_empty_idx,class_idx),2));
        %empty_query_cluster_1
%         temp_QueryFea = getqueryfeature(temp_QueryFea, sort_cluster_idx(class_idx), type,alpha,t_dist,normweight);
        
        ClickCount = sum(sort_feature_click(query_g1_empty_idx,class_idx),2);
        temp_QueryFea = getqueryfeature(sort_feature_click(query_g1_empty_idx,class_idx), sort_cluster_idx(class_idx), type,alpha,t_dist);
        [~, cluster_idx_empty] = emptyQueryClustering(temp_QueryFea, query_g1_empty_idx,Keachone_query,len_empty,gamma,P, normweight);
                
        
        cluster_idx_empty = cluster_idx_empty + max(KNNIDX);
        KNNIDX(query_g1_empty_idx) = cluster_idx_empty;
    elseif type_empty_cluster == 2
        Neachone_query_empty = cellfun(@(x) numel(x),each_class_query_g1_empty_idx);
        Keachone_query_empty = getKeachclass_query(method,k_cluster_Q,Neachone_query_empty);
        for i = 1:len_empty
            h = tic;
            if ~isempty(each_class_query_g1_empty_idx{i})
                Q_index_E = each_class_query_g1_empty_idx{i};
  
                range = 1:length(Q_index_E);
                if type == 'WP'
                    t_dist = EuDist2(data_fea(data_label==each_class_index_query_g1_empty(i),:),data_fea(data_label==each_class_index_query_g1_empty(i),:));
                else
                    t_dist = [];
                end
                %%%pro
                each_class_idx = sum(Neachone_img(1:each_class_index_query_g1_empty(i)))-Neachone_img(each_class_index_query_g1_empty(i))+1:sum(Neachone_img(1:each_class_index_query_g1_empty(i)));
                
%                 %%%normal
%                 QueryFea = sort_feature_click(Q_index_E,each_class_idx) ./ repmat(sqrt(sum(sort_feature_click(Q_index_E,each_class_idx).^2,2)),1,size(sort_feature_click(Q_index_E,each_class_idx),2));
%                 QueryFea = getqueryfeature(QueryFea, sort_cluster_idx(each_class_idx), type,alpha,t_dist,normweight);

                ClickCount = sum(sort_feature_click(Q_index_E,each_class_idx),2);
                QueryFea = getqueryfeature(sort_feature_click(Q_index_E,each_class_idx), sort_cluster_idx(each_class_idx), type,alpha,t_dist);
                
                %%clustering
                temp_QueryFea = QueryFea';
                cluster_data = full(temp_QueryFea);
                cluster_number = Keachone_query_empty(i);
                [cluster_vector ,cluster_idx_q] = Myclustering(cluster_data,...
                    cluster_number,i,cluster_method,gamma,P,theta,T, textweight,normweight,ClickCount);
                
                
       
                cluster_idx_q = cluster_idx_q + max(KNNIDX);

                %%1nn
                [dist, idx] = min(EuDist2(QueryFea(range,:), cluster_vector'), [], 2);
                idx = double(cluster_idx_q(idx'))';
                tmp= KNNDis(Q_index_E) < dist;
                KNNIDX(Q_index_E) =  (tmp) .* KNNIDX(Q_index_E) + ~tmp .* idx;
        %         KNNDis(range2) =  (tmp) .* KNNDis(range2) + ~tmp .* dist;
                KNNDis(Q_index_E(tmp==0)) = dist(tmp==0);
    %             IDX(Q_index_E) = KNNIDX;
                disp(i);
            end
            h1 = toc(h);
            sprintf('%d : %d',i,h1);
        end
    end
    
    h1= toc(h);
    sprintf(' emptyQueryClustering %d', h1)
    whos
    
    IDX = KNNIDX;
    
%     eval(['save(''' sname ''', ''IDX'')']);
    save(fullfile('result', [sname, 'c-', num2str(type_empty_cluster), '-IDX.mat']), 'IDX');
    disp(['save successful:' sname]);
end
