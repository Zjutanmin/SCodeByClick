%根据测试query的k个近邻query的clickcount来求权值

path = pwd;
cd ..;

load(fullfile(pwd, 'VisualSimilarity/SampleClick.mat'));
load(fullfile(pwd, 'VisualSimilarity/New_SampleImageList.mat'));
load(fullfile(pwd, 'VisualSimilarity/UniqueSampleImage.mat'));
load(fullfile(pwd, 'VisualSimilarity/SampleQuery.mat'));

cd(path);

len = length(UniqueSampleImage);
clickvector = cell(length(New_SampleImageList),1);

for i = 1:length(New_SampleImageList)
    v = zeros(len,1);
    v(New_SampleImageList{i}) = SampleClick{i};
    clickvector{i} = sparse(v);
    disp(i);
end

save('click_vector.mat','clickvector');

computedis_clickvector;



