

% len = length(SampleQuery);
len = 123006;
% similar_clickvector = cell(len,1);
similar_clickvector = zeros(len,len);
tic;
for i = 1:len
    for j = i+1:len
        similar_clickvector(i,j) = norm(clickvector{i}-clickvector{j});
    end
    disp(i);
end
toc;
similar_clickvector = sparse(similar_clickvector);

save('similar_clickvector.mat','similar_clickvector');






