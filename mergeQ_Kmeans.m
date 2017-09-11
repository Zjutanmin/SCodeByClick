function IDX = mergeQ_Kmeans(theta)

addpath(genpath('home/weichen/'))
addpath(genpath('/home/weichen/computeWeight/ksvd+omp/ksvdbox'))
addpath(genpath('/home/weichen/computeWeight/ksvd+omp/ompbox'))
addpath('/home/weichen/computeWeight/1-283')
method = 'N';
k_cluster_I = 20;
k_cluster_Q = 20;
type = 'WP';
alpha = 0.4;
P = 5;
hotrate = 0.4;  T = 0.04;

sname = [method,'-',type,'-k-', num2str(k_cluster_I),'-',num2str(k_cluster_Q),'-alpha-',num2str(alpha),'-theta-',num2str(theta),'-T-',num2str(T),'-P-',num2str(P),'.mat'];
addpath(genpath('vlfeat'))
addpath('/home/Weichen/acc')

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

    
    h1= toc(h);
    sprintf('load %d', h1)
    whos
    h = tic;
    
    try
        load(['kImgC/' 'ImgC-', method, '-', num2str(k_cluster_I), '.mat'], 'cluster_idx')
    catch
        Keachone_img = getKeachclass_img(method, k_cluster_I, Neachone_img);
        cluster_idx = ImageClustering(Keachone_img);
        save(['kImgC/' 'ImgC-', method, '-', num2str(k_cluster_I), '.mat'], 'cluster_idx')
        disp('save successful');
    end
    
    h1= toc(h);
    sprintf('ImagClustering %d', h1)
    whos
    h = tic;
    

    [sort_cluster_idx, sort_idx] = sort(cluster_idx);
    sort_feature_click = feature_click(:,sort_idx);

      
    Keachone_query = getKeachclass_query(method, k_cluster_Q, Neachone_query);
    
    if type == 'WP'
        load('image_click_Dog283_0_CNN_Alex1_ND_S_S1_data_normal.mat');

    end
    
    
    try 
        load(['kKNN/' 'KNN-',method,'-',num2str(k_cluster_I),'-',num2str(k_cluster_Q),'-',type,'-',num2str(alpha),'-',num2str(P),'.mat'],'KNNIDX','KNNDis');
    catch    
    KNNIDX = ones(length(each_query_class_index), 1) * 0;
    KNNDis = ones(length(each_query_class_index), 1) * inf;

    for i = 1 : length( Keachone_query)
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
          
            each_class_idx = sum(Neachone_img(1:i))-Neachone_img(i)+1:sum(Neachone_img(1:i));
         
            QueryFea = sort_feature_click(Q_index_E,each_class_idx);
            
            tic
            QueryFea = getqueryfeature(QueryFea, sort_cluster_idx(each_class_idx), type,alpha,t_dist);
            toc
            
            %%clustering
            temp_QueryFea = QueryFea(range1,:)';
%             [len_clickvector,~ ] = size(temp_QueryFea);
       
       matx = zeros(size(temp_QueryFea,2));
        for ii = 1 : size(temp_QueryFea,2)
        sumx = sum(temp_QueryFea(:,ii));
        if(sumx == 0)
            sumx = 1;
        end
        
        matx(ii,ii) = 1/sumx;
        end
        if(find(temp_QueryFea))
            norm_temp_QueryFea = temp_QueryFea * matx;
        end
            
% %             if(find(isnan(temp_QueryFea)))
% %                 temp_QueryFea(find(isnan(temp_QueryFea))) = 0;
% %             end
% %              if(find(isnan(norm_temp_QueryFea)))
% %                 norm_temp_QueryFea(find(isnan(norm_temp_QueryFea))) = 0;
% %             end
            Nc = floor(size(temp_QueryFea,2)*hotrate/P);
            if(Nc==0)
                Nc = size(temp_QueryFea,2);
            end
           ts_name = [num2str(i),'.mat'];
           load(ts_name);
           if(~isempty(clus))
            [~,kmIDX] = kmeans(full(clus),Nc);
           else
           [~,kmIDX] = vl_kmeans(norm_temp_QueryFea,Nc);
           end
          tempcount = [];D = [];
          for ai = 1 : size(temp_QueryFea,2)
              tempcount(1,ai) = sum(temp_QueryFea(:,ai));
          end
          
       
          params.initdict = D;
         params.dictsize = size(D,2);
         
         
        pX = full(norm_temp_QueryFea);      
       
        
        addpath(genpath('/home/Weichen/code/ksvd+omp'));       
        
       
       
        params.data = pX;           % %��ݱ���?
        params.Tdata = size(pX,2)*hotrate*T;           % %Լ�������Լ��ֵ
        params.fixedDic = 1;
        params.iternum = 30;        % % ������
        params.memusage = 'high';
        params.fixedcode = 0;          % %���벻��
        params.method_dic = 'L0';      % %����0���Ƿ���1�����ַ�ʽ
        params.SCode.Code = [];  %%%Min Tan
        params.lamda = 1;  % %��Ҫ���ڣ�����Ϊ1ʱ��Ĳ�?
        
 
         
          if(size(temp_QueryFea,2)>theta)   
             [cluster_vector,cluster_idx_q ] = vl_kmeans(pX,Nc);  
           
        
         elseif(size(temp_QueryFea,2)>Nc)
             [~,IDX ] = vl_kmeans(pX,Nc); 
            [a,b,c] = unique(IDX);
            pD = [];
            code_size = [];code_size1 = [];code_zise2 = [];
            
             for pi = 1:length(a)
                 index = find(c == pi);
                 [~, iind] = sort(tempcount( index), 'descend');
                 %iind = iind(1:max(1, floor(length(index)*hotrate)));                 
                 pD = [pD,pX(:, index(iind)) ];
                 code_size(pi) = length(iind);
                 code_size1(pi) = size(pD,2);
                 
             end

	    code_size2(1) = 0;
             
             for pi = 2:length(code_size1)+1
             code_size2(pi) = code_size1(pi-1);
             end

             params.Tdata = T;
             params.dictsize = size(pD,2);
             params.initdict = pD;
                 [~ ,code] = ksvd(params,'');
                 code1 = [];
                 
              for pi = 1:length(a)
                 start1 = 1+code_size2(pi);
                  end1 = code_size2(pi+1);
                  code1 = [code1;sum(      code(start1:end1,:)  ,1)]  ;
              end
             [~,cluster_idx_q] = max(code1, [], 1);     
              
              cluster_vector = params.data(:,code_size1);
                 
           
             
         else  
            code = eye(size(pX,2));
            cluster_vector = params.initdict;
             [~,cluster_idx_q] = max(code,[],1);
         end
    
% %        pathname = ['/home/weichen/DICT'];
% %        
% %        filename = ['/Dictionary-',method,'-',num2str(k_cluster_I),'-',num2str(k_cluster_Q),'-',type,'-',num2str(alpha),'-P-',num2str(P),'-T-',num2str(theta),'-dpp-',num2str(isdpp),'.mat'];   
% %        pfname = [pathname,filename];
% %        save(pfname,'cluster_vector');
            
            cluster_idx_q = cluster_idx_q + max(KNNIDX);
    
            KNNIDX(Q_index_E(range1)) = cluster_idx_q;

            
%             cluster_vector = cluster_vector(len_clickvector,:);


            %%1nn
            range = Q_index_E(range2);
            [dist, idx] = min(EuDist2(QueryFea(range2,:), cluster_vector'), [], 2);
            idx = double(cluster_idx_q(idx'))';
            tmp= KNNDis(range) < dist;
            KNNIDX(range) =  (tmp) .* KNNIDX(range) + ~tmp .* idx;
  
            KNNDis(range(tmp==0)) = dist(tmp==0);

            disp(i);
        end
        h1 = toc(h);
        sprintf('%d : %d',i,h1);
    end
    h1= toc(h);
    sprintf(' QueryClustering %d', h1)
    whos
    h = tic;
    
    save(['kKNN/' 'KNN-',method,'-',num2str(k_cluster_I),'-',num2str(k_cluster_Q),'-',type,'-',num2str(alpha),'-',num2str(P),'-',num2str(theta),'.mat'],'KNNIDX','KNNDis');
    disp('save successful');
    end    
    
   
    
    h1= toc(h);
    sprintf(' emptyQueryClustering %d', h1)
    whos
           
    IDX = KNNIDX;
    

    save(['kresult/' sname],'IDX');
    disp(['save successful:' sname]);
