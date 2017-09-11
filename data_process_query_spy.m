
clear all;
load('each_query_class.mat');
% load('each_class_query.mat')
load('query123_new.mat');

% idx_1 = unique(cell2mat(Q_index_E1));
% idx_2 = unique(cell2mat(Q_index_E2'));
% idx_3 = unique(cell2mat(Q_index_E3));

% each_query = each_query_class_index;
% each_query123(idx_1,1) = 1;
% each_query123(idx_2,1) = 2;
% each_query123(idx_3,1) = 3;

randquery_idx_1 = cell(length(Q_index_E1),1);
randquery_idx_2 = cell(length(Q_index_E2),1);
randquery_idx_3 = cell(length(Q_index_E3),1);

num = 5;

for i = 1:length(Q_index_E1)
    if ~isempty(Q_index_E1{i})
        s = numel(Q_index_E1{i});
        if s>num
            i_idx = randperm(s,num);
            randquery_idx_1{i} = Q_index_E1{i}(i_idx);
        else
            randquery_idx_1{i} = Q_index_E1{i};
        end
    end
end

for i = 1:length(Q_index_E2)
    if ~isempty(Q_index_E2{i})
        s = numel(Q_index_E2{i});
        if s>num
            i_idx = randperm(s,num);
            randquery_idx_2{i} = Q_index_E2{i}(i_idx)';
        else
            randquery_idx_2{i} = Q_index_E2{i}';
        end
    end
end

for i = 1:length(Q_index_E3)
    if ~isempty(Q_index_E3{i})
        s = numel(Q_index_E3{i});
        if s>num
            i_idx = randperm(s,num);
            randquery_idx_3{i} = Q_index_E3{i}(i_idx);
        else
            randquery_idx_3{i} = Q_index_E3{i};
        end    
    end
end

randquery_idx_1 = cell2mat(randquery_idx_1);
randquery_idx_2 = cell2mat(randquery_idx_2);
randquery_idx_3 = cell2mat(randquery_idx_3);

each_query_idx = [randquery_idx_1;randquery_idx_2;randquery_idx_3];


len1 = length(randquery_idx_1);
len2 = length(randquery_idx_2);
len3 = length(randquery_idx_3);

each_query = cell(len1+len2+len3,2);
each_query(1:len1,1) = each_query_class_index(randquery_idx_1);
each_query(1+len1:len1+len2,1) = each_query_class_index(randquery_idx_2);
each_query(1+len1+len2:len1+len2+len3,1) = each_query_class_index(randquery_idx_3);


each_query(1:len1,2) = {1};
each_query(1+len1:len1+len2,2) = {2};
each_query(1+len1+len2:len1+len2+len3,2) = {3};


