function [cluster_vector,cluster_idx_q] = Myclustering(temp_QueryFea,...
    Keachone_query,i,cluster_method,gamma,P,theta,T, textweight,normweight,ClickCount)

temp_QueryFea = Mynomalize(temp_QueryFea, 1, normweight);

global setting;

%%KMENA
if setting.UsePGNc == 1  %all for Keachone_query(i)
    Nc1 = floor(size(temp_QueryFea,2)*gamma/P);
    if(Nc1==0)
        Nc1 = size(temp_QueryFea,2);
    end
    Nc2 = Nc1;
else
    Nc1 = Keachone_query(i);
    Nc2 = Nc1;
end



switch abs(cluster_method)
    case 1  %kmeans
        [cluster_vector,cluster_idx_q] = vl_kmeans(full(temp_QueryFea),Nc1);
    case 2  %Hot
        [cluster_vector,cluster_idx_q] =Sparse_cluster(full(temp_QueryFea),i,gamma,P,theta,T, 1, textweight, sign(cluster_method),Nc2,normweight,ClickCount);
    case 3  %Ksvd
        [cluster_vector,cluster_idx_q] =Sparse_cluster(full(temp_QueryFea),i,gamma,P,theta,T, 0, textweight, sign(cluster_method),Nc2,normweight,ClickCount);
end
        
        
function   [cluster_vector,cluster_idx_q] = Sparse_cluster(temp_QueryFea,i,...
        gamma,P,theta,T,type,textweight, Code_method,Nc,normweight,tempcount)
if(i~=116)
ts_name = [num2str(i),'.mat'];
load(ts_name);
end
if(i==116)
    load('116_1.mat');
    load('116_2.mat');
    load('116_3.mat');
    clus = [clus1,clus2,clus3];
end
clus = bsxfun(@times, full(clus), 1./ sum(full(clus),1));
% % if(~isempty(clus))
% %     [~,kmIDX] = vl_kmeans(full(clus),Nc); %%%text feature by similarity
% % else
% %     [~,kmIDX] = vl_kmeans(norm_temp_QueryFea,Nc); %%%click feature by similarity
% % end

if textweight == 0
else if textweight == 1
        QueryFea = textweight * clus;
    else
        if normweight == 1
            QueryFea = [textweight * clus; (1-textweight) * temp_QueryFea];
        else if normweight == 2
                QueryFea = [sqrt(textweight) * clus; sqrt(1-textweight) * temp_QueryFea];
            end
        end
    end
end

[~,kmIDX] = vl_kmeans(QueryFea,Nc);%%%text feature by similarity

% tempcount = [];
D = [];
% for ai = 1 : size(temp_QueryFea,2)
%     tempcount(1,ai) = sum(temp_QueryFea(:,ai));
% end
%search for Hot dictionary
D = [];D_Label = [];
for kmi = 1 : Nc
    kindex = find(kmIDX == kmi);
    [~,kia] = sort(tempcount(kindex),'descend');
    pk = min(P,length(kia));kia = kia(1:pk);
    D = [D,temp_QueryFea(:,kindex(kia))];
   % tempcount(kindex(kia))
    D_Label = [D_Label; ones(pk,1)*kmi];
end


pX = full(temp_QueryFea);
params.initdict = normcols(D);
if(size(temp_QueryFea,2)>theta)  %%%for sparse coding 
    addpath(genpath('ksvd+omp'))
    
%     params.initdict = D;
    params.dictsize = size(D,2);
    addpath(genpath('/home/Weichen/code/ksvd+omp'));   
    params.data = pX;  
    params.Tdata = ceil(size(D,2)*T);
    params.fixedDic = type;
    params.iternum = 20;     
    params.memusage = 'high';
    params.fixedcode = 0;    
    params.method_dic = 'L0';  
    params.SCode.Code = [];  %%%Min Tan
    params.lamda = 1; 
    
    [D_New,code] = ksvd(params,'');code = full(code);
    
    
    if (type == 0)  && (P> 1 )%%%%KSVD, then update D_Label
        [~,D_Label] = vl_kmeans(D_New,Nc); %%%text feature by similarity
    end
    
    
    
    Nc_cell = mat2cell([1 : Nc]', ones(Nc,1), 1);
    if Code_method == 1
        code = cell2mat(cellfun(@(x) sum(code(find(D_Label==x),:),1), Nc_cell,'UniformOutput',0));
    else
        code = cell2mat(cellfun(@(x) mean(code(find(D_Label==x),:),1), Nc_cell,'UniformOutput',0));
    end
    
    cluster_vector = D_New;
    [~,cluster_idx_q] = max(code,[],1);

else
    if(size(temp_QueryFea,2)>Nc)
        [cluster_vector,cluster_idx_q ] = vl_kmeans(pX,Nc);  
    else
        code = eye(size(pX,2));
        cluster_vector = params.initdict;
        [~,cluster_idx_q] = max(code,[],1);
    end
end
end

end    
% function   [cluster_vector,cluster_idx_q] = Hot_cluster(temp_QueryFea,i,gamma,P,theta,T)
%      matx = zeros(size(temp_QueryFea,2));
%         for ii = 1 : size(temp_QueryFea,2)
%         sumx = sum(temp_QueryFea(:,ii));
%         if(sumx == 0)
%             sumx = 1;
%         end
%         
%         matx(ii,ii) = 1/sumx;
%         end
%         if(find(temp_QueryFea))
%             norm_temp_QueryFea = temp_QueryFea * matx;
%         end
%             
%             if(find(isnan(temp_QueryFea)))
%                 temp_QueryFea(find(isnan(temp_QueryFea))) = 0;
%             end
%              if(find(isnan(norm_temp_QueryFea)))
%                 norm_temp_QueryFea(find(isnan(norm_temp_QueryFea))) = 0;
%             end
%             Nc = floor(size(temp_QueryFea,2)*gamma/P);
%             if(Nc==0)
%                 Nc = size(temp_QueryFea,2);
%             end
%            ts_name = [num2str(i),'.mat'];
%            load(ts_name);
%            if(~isempty(clus))
%             [~,kmIDX] = vl_kmeans(full(clus),Nc);
%            else
%            [~,kmIDX] = vl_kmeans(norm_temp_QueryFea,Nc);
%            end
%           tempcount = [];D = [];
%           for ai = 1 : size(temp_QueryFea,2)
%               tempcount(1,ai) = sum(temp_QueryFea(:,ai));
%           end
%          
%           
%           %search for Hot dictionary
%            for kmi = 1 : Nc
%                kindex = find(kmIDX == kmi);
%                [~,kia] = sort(tempcount(kindex),'descend');
%                D = [D,temp_QueryFea(:,kia(1:min(P,length(kia))))];
%            end
%        
%         
%           
%           params.initdict = D;          
%          params.dictsize = size(D,2);
%          
%          
%         pX = full(norm_temp_QueryFea);
%        
%         
%         addpath(genpath('/home/Weichen/code/ksvd+omp'));       
%         
%        
%        
%         params.data = pX;          
%         params.Tdata = size(pX,2)*gamma*T;       
%         params.fixedDic = 1;
%         params.iternum = 30;     
%         params.memusage = 'high';
%         params.fixedcode = 0;       
%         params.method_dic = 'L0';     
%         params.SCode.Code = [];  %%%Min Tan
%         params.lamda = 1; 
%         
%          
%       if(size(temp_QueryFea,2)>theta)   
%                 
%             [~,IDX ] = vl_kmeans(pX,Nc); 
%             [a,~,c] = unique(IDX);
%             pD = [];
%             code_size = [];code_size1 = [];code_zise2 = [];
%             
%              for pi = 1:length(a)
%                  index = find(c == pi);
%                  [~, iind] = sort(tempcount( index), 'descend');              
%                  pD = [pD,pX(:, index(iind)) ];
%                  code_size(pi) = length(iind);
%                  code_size1(pi) = size(pD,2);
%                  
%              end
% 
% 	    code_size2(1) = 0;
%              
%              for pi = 2:length(code_size1)+1
%              code_size2(pi) = code_size1(pi-1);
%              end
% 
%              params.Tdata = T;
%              params.dictsize = size(pD,2);
%              params.initdict = pD;
%               [~ ,code] = ksvd(params,'');
%              code1 = [];
%                  
%               for pi = 1:length(a)
%                  start1 = 1+code_size2(pi);
%                   end1 = code_size2(pi+1);
%                   code1 = [code1;sum(      code(start1:end1,:)  ,1)]  ;
%               end
%              [~,cluster_idx_q] = max(code1, [], 1);     
%               
%               cluster_vector = params.data(:,code_size1);
%                  
%            
%            
%         
%          elseif(size(temp_QueryFea,2)>Nc)
%               [cluster_vector,cluster_idx_q ] = vl_kmeans(pX,Nc);  
%          else  
%             code = eye(size(pX,2));
%             cluster_vector = params.initdict;
%              [~,cluster_idx_q] = max(code,[],1);
%     end        
% end
% 
% function   [cluster_vector,cluster_idx_q] = KH_cluster(temp_QueryFea,i,gamma,P,theta,T)
%      matx = zeros(size(temp_QueryFea,2));
%         for ii = 1 : size(temp_QueryFea,2)
%         sumx = sum(temp_QueryFea(:,ii));
%         if(sumx == 0)
%             sumx = 1;
%         end
%         
%         matx(ii,ii) = 1/sumx;
%         end
%         if(find(temp_QueryFea))
%             norm_temp_QueryFea = temp_QueryFea * matx;
%         end
%             
%             if(find(isnan(temp_QueryFea)))
%                 temp_QueryFea(find(isnan(temp_QueryFea))) = 0;
%             end
%              if(find(isnan(norm_temp_QueryFea)))
%                 norm_temp_QueryFea(find(isnan(norm_temp_QueryFea))) = 0;
%             end
%             Nc = floor(size(temp_QueryFea,2)*gamma/P);
%             if(Nc==0)
%                 Nc = size(temp_QueryFea,2);
%             end
%            ts_name = [num2str(i),'.mat'];
%            load(ts_name);
%            if(~isempty(clus))
%             [~,kmIDX] = vl_kmeans(full(clus),Nc);
%            else
%            [~,kmIDX] = vl_kmeans(norm_temp_QueryFea,Nc);
%            end
%           tempcount = [];D = [];
%           for ai = 1 : size(temp_QueryFea,2)
%               tempcount(1,ai) = sum(temp_QueryFea(:,ai));
%           end
%          
%           
%           %search for Hot dictionary
%            for kmi = 1 : Nc
%                kindex = find(kmIDX == kmi);
%                [~,kia] = sort(tempcount(kindex),'descend');
%                D = [D,temp_QueryFea(:,kia(1:min(P,length(kia))))];
%            end
%        
%         
%           
%           params.initdict = D;          
%          params.dictsize = size(D,2);
%          
%          
%         pX = full(norm_temp_QueryFea);
%        
%         
%         addpath(genpath('/home/Weichen/code/ksvd+omp'));       
%         
%        
%        
%         params.data = pX;          
%         params.Tdata = size(pX,2)*gamma*T;       
%         params.fixedDic = 0;
%         params.iternum = 30;     
%         params.memusage = 'high';
%         params.fixedcode = 0;       
%         params.method_dic = 'L0';     
%         params.SCode.Code = [];  %%%Min Tan
%         params.lamda = 1; 
%         
%          
%           if(size(temp_QueryFea,2)>theta)   
%                 
%             [~,IDX ] = vl_kmeans(pX,Nc); 
%             [a,~,c] = unique(IDX);
%             pD = [];
%             code_size = [];code_size1 = [];code_zise2 = [];
%             
%              for pi = 1:length(a)
%                  index = find(c == pi);
%                  [~, iind] = sort(tempcount( index), 'descend');              
%                  pD = [pD,pX(:, index(iind)) ];
%                  code_size(pi) = length(iind);
%                  code_size1(pi) = size(pD,2);
%                  
%              end
% 
% 	    code_size2(1) = 0;
%              
%              for pi = 2:length(code_size1)+1
%              code_size2(pi) = code_size1(pi-1);
%              end
% 
%              params.Tdata = T;
%              params.dictsize = size(pD,2);
%              params.initdict = pD;
%              [~ ,code] = ksvd(params,'');
%              code1 = [];
%                  
%               for pi = 1:length(a)
%                  start1 = 1+code_size2(pi);
%                   end1 = code_size2(pi+1);
%                   code1 = [code1;sum(      code(start1:end1,:)  ,1)]  ;
%               end
%              [~,cluster_idx_q] = max(code1, [], 1);     
%               
%               cluster_vector = params.data(:,code_size1);
%                  
%            
%            
%         
%          elseif(size(temp_QueryFea,2)>Nc)
%               [cluster_vector,cluster_idx_q ] = vl_kmeans(pX,Nc);  
%          else  
%             code = eye(size(pX,2));
%             cluster_vector = params.initdict;
%              [~,cluster_idx_q] = max(code,[],1);
%     end        
% end
%  
% % % % 
% % % % function [C,D] = vl_kmeans(A,B) %%%text feature by similarity
% % % % C = A(:, 1:B);
% % % % D = [[1:B], ceil(rand(1,size(A,2) - B)*B)];
% % % % end
% % % % 
% % % % function [C,D] = kmeans(A,B) %%%text feature by similarity
% % % % C = A(:, 1:B);
% % % % D = [[1:B], ceil(rand(1,size(A,2) - B)*B)];
% % % % end
