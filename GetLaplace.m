function Lap = GetLaplace(data_N, LapName, type, K)
if nargin < 3
    type = 'N';
    K = 20;
end
Lapstr = [type '-' num2str(K) '-s.mat'];
ImageN = [LapName '-ImgC-' Lapstr];

try
    load(fullfile('SimMatrix', ImageN), 'G');
catch
    G = GetSimMatrix(data_N, type, K);
    save(fullfile('SimMatrix', ImageN), 'G');
end
D = diag(sum(G, 2));
Lap = D-G;





function each_class_imageCNN_s = GetSimMatrix(ddata_N, type, k)
% method = {'N','P'};
% c_k_cluster_I = {[10, 20, 50, 100], [5, 10, 50, 100, 500, 1000]};
% 
% Mrange = 1:length(method);

% load('ddata_N.mat');
% data_fea_N = myNormlize(data_fea,1);

% if ~exist('Qrange', 'var') 
%     Qrange = 1:length(c_k_cluster_I_T{Mrange{i}});
% end

% for i_method = Mrange
%     k_cluster_I = c_k_cluster_I{i_method};
%     Qrange = 1:length(k_cluster_I);
%     for i_k_cluster_I = Qrange;
        tic;
%         method{i_method}
        
%         sname = ['ImgC/ImgC-' method{i_method} '-' num2str(k_cluster_I(i_k_cluster_I)) '.mat'];
        sname = ['ImgC/ImgC-' type '-' num2str(k) '.mat'];
        
        load(sname);
        
        each_class_image_idx = arrayfun(@(x) find(cluster_idx == x),1:max(cluster_idx),'UniformOutput',false);
        each_class_image_vector = cellfun(@(x) ddata_N(x,:),each_class_image_idx,'UniformOutput',false);
% % %         each_class_imageCNN_d = cellfun(@(x) EuDist2(x,x),each_class_image_vector,'UniformOutput',false);
        each_class_imageCNN_d = cellfun(@getSim,each_class_image_vector,'UniformOutput',false);
        
        
        x_idx = cell(length(each_class_imageCNN_d),1);
        y_idx = cell(length(each_class_imageCNN_d),1);
        s = cell(length(each_class_imageCNN_d),1);
        for i = 1:length(each_class_imageCNN_d)
            len = length(each_class_image_idx{i});
            x_idx{i,1} = repmat(each_class_image_idx{i},len,1);
            y_idx{i,1} = repmat(each_class_image_idx{i},1,len)';
            y_idx{i,1} = y_idx{i,1}(:);
            s{i,1} = each_class_imageCNN_d{i}(:);
        end
        x_idx = cell2mat(x_idx);
        y_idx = cell2mat(y_idx);
        s = cell2mat(s);
        each_class_imageCNN_s = sparse(x_idx,y_idx,s,length(cluster_idx),length(cluster_idx));
        clear x_idx;
        clear y_idx;
        clear s;
% % % %         [x_idx,y_idx,s] = find(each_class_imageCNN_d);
% % % %         each_class_imageCNN_s = sparse(x_idx,y_idx,s,length(cluster_idx),length(cluster_idx));
        
%         each_class_imageCNN_s = zeros(length)
% % %         for i = 1:length(each_class_imageCNN_d)
% % %             each_class_imageCNN_s(each_class_image_idx{i},each_class_image_idx{i}) = each_class_imageCNN_d{i};
% % %         end
%         save(['ImgC/ImgC-' type '-' num2str(k) '-s.mat'],'each_class_imageCNN_s')
        disp(['save ' type '-' num2str(k)]);
        toc;
%     end
% end
