%根据测试query的clickcount来求权值

path = pwd;
cd ..;

load(fullfile(pwd, 'VisualSimilarity/SampleClick.mat'));

weight_y = cell(length(SampleClick),1);

for i = 1:length(SampleClick)
    s = sum(SampleClick{i});
    weight_y{i} = SampleClick{i}*length(SampleClick{i})/s;
end

cd(path);
save('weight_y.mat','weight_y');
