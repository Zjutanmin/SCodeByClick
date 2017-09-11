function getclickfeature1(data_dir, fea_dir, database, featparaclick, featpara)
if nargin < 4
    featparaclick = '';
    featpara = '';
end


        
idx = findstr(data_dir, featparaclick);


try
    load([data_dir(1:idx(1)-1), data_dir(idx(1)+length(featparaclick):end)], 'data_fea', 'data_label')
catch




try
load([data_dir(1:idx(1)-1), data_dir(idx(1)+length(featparaclick):end)], 'data_fea', 'data_label')

% load(fullfile(fea_dir, 'Dog_ClickImg'), 'img_Code', 'img_ID', 'img_name')

catch
    
fea_dir = 'D:\m\data\dog';
load(fullfile(fea_dir, ['DogData.mat']));


cellinfo = {};
imgfidinfo = {};
imgtidinfo = {};
NewName = {};

temp = 0;

for i = 1:length(Dog_Breed)
    for j = 1:length(Dog_Image{i})
        temp = temp +1 ;
        cellinfo{end+1} = Dog_Breed{i};
        imgfidinfo{end+1} = j;
        imgtidinfo{end+1} = temp;
        NewName{end+1} = fullfile(Dog_Breed{i}, [num2str(j),'.jpg']);
    end
    inum(i) = length(Dog_Image{i});
end

index = cellfun(@(x) findcstr(x, Dog_Breed), database.cname');
inum = [0, cumsum(inum)]; inum = inum(1:end-1);
sindex = inum(index);sindex = sindex';
[~,B,~]=cellfun(@fileparts, database.orgimpath, 'UniformOutput', false);
B= cellfun(@str2num, B, 'UniformOutput', false);
B = cell2mat(B');
fid = sindex(database.label)+B;

load(fullfile(fea_dir, 'Sort_Dog_ID_sindex.mat'), 'ic')
load(fullfile(fea_dir, 'US_DogClick_Q.mat'), 'Dog_QQ');
load(fullfile(fea_dir, 'sparse.mat'), 'D')
load(fullfile(fea_dir, 'Dog_ClickFG'), 'Dog_ClickFG')

Fea_N = sparse(temp, length(Dog_QQ));
% Fea = sparse(temp, length(Dog_QQ));
Fea_T = zeros(temp, 1);
Fea_N(ic,:) = D;
Fea_T(ic) = cell2mat(Dog_ClickFG(:, 2));
% Fea(ic,:) = bsxfun(@times, D, cell2mat(Dog_ClickFG(:, 2)));

load(fullfile(fea_dir, 'Sort_Dog_ID.mat'), 'Sort_Dog_ID')
img_ID = cell(temp, 1);
img_ID(ic,:) = Sort_Dog_ID;

load(fullfile(fea_dir, 'Sort_Dog_Code.mat'), 'Sort_Dog_Code')
img_Code = cell(temp, 1);
img_Code(ic,:) = Sort_Dog_Code;

Fea_N = Fea_N(fid,:);
Fea_T = Fea_T(fid,:);
% Fea = Fea(fid,:);
img_ID = img_ID(fid,:);
img_Code = img_Code(fid,:);

index = find(sum(Fea_N, 1));

Dog_QQ_N = Dog_QQ(index);
img_Fea_N = Fea_N(:, index);
img_Fea_T = Fea_T;
% img_Fea = Fea(:, index);
img_label = (database.label);
img_name = (database.orgimpath);
save(fullfile(fea_dir(1:idx(1)-1), 'Dog_ClickFea'), 'Dog_QQ_N', 'img_Fea_N', ...
    'img_Fea_T', 'img_label', 'img_name')
save(fullfile(fea_dir(1:idx(1)-1), 'Dog_ClickImg'), 'img_Code', 'img_ID', 'img_name')
end


try
    load(fea_dir(1:idx(1)-1), 'data_label')
catch
    data_fea = img_Fea_N;
    data_label = img_label;
    save(fea_dir(1:idx(1)-1), 'data_fea', 'data_label')
end

save([data_dir(1:idx(1)-1), data_dir(idx(1)+length(featparaclick):end)], 'data_fea', 'data_label')
end

[N, M] = size(data_fea);
try
    load(data_dir, 'data_fea', 'data_label')
catch
    method = featpara{2};
    k_cluster_I = featpara{4};
    k_cluster_Q = featpara{5};
    type = featpara{3};
    alpha = featpara{6};
     IDX = mergeQ2(method, k_cluster_I, k_cluster_Q, type, alpha);
     NClass = length(unique(IDX));
     data_fea_T = sparse(size(data_fea, 1), NClass);
     for j = NClass
         data_fea_T(:, j) = sum(data_fea(:, find(IDX == j)), 2);
     end
    data_fea = data_fea_T;
    save(data_dir, 'data_fea', 'data_label')

end
