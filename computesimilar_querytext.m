

origin_path = pwd;
cd ..
path = fullfile(pwd, 'TextSimilarity');

cd(path);
load('StopWords.mat');
% len = length(SampleQuery);
len =200;
similar_querytext = zeros(len,len);

tic
for i = 1:length(SampleQuery)
    for j = 1:length(SampleQuery)
        similar_querytext(i,j) = Fun_TextSimilarity(SampleQuery{i},SampleQuery{j});
        disp(j);
    end
    disp(i);
end
toc

save('similar_querytext.mat','similar_querytext');

