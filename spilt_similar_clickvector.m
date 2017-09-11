

for i = 20000:20000:120000
    tic;
    eval(['similar_clickvector' num2str(i) '= similar_clickvector(' num2str(i-19999) ':' num2str(i) ',:);']);
    eval(['save(''similar_clickvector' num2str(i) '.mat'', ''similar_clickvector' num2str(i) ''');']);
    eval(['clear similar_clickvector' num2str(i) ';']);
    disp(i);
toc;
end