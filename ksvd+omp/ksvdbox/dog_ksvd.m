addpath('D:\m\code\ksvd+omp\ksvdbox')
addpath('D:\m\code\ksvd+omp\ompbox')
load('Dict.mat');
load('image_click_Dog283_0_img_Fea_Clickcount.mat')

% % X = img_Fea_Clickcount';

X = img_Fea_Clickcount';

k = 4;
m = size(D,1);
params.lamdaRange = [0.1:0.1:1];

params.data = X;           % %数据本身
params.Tdata = k;           % %约束项范数后的约束值
params.dictsize = m;        % %字典需要基向量的个数(e1,e2.....em)
params.iternum = 30;        % % 迭代次数
params.memusage = 'high';
params.fixedcode = 0;          % %编码不变
params.method_dic = 'L0';      % %范数0还是范数1的两种方式
params.SCode.Code = [];  %%%Min Tan
params.lamda = 1;         % %需要调节，范数为1时候的参数

params.fixedDic =0;
[~,code] = ksvd(params,'');    %%%%入口

%% show results %%

figure; plot(err); title('K-SVD error convergence');
xlabel('Iteration'); ylabel('RMSE');

printf('  Dictionary size: %d x %d', n, m);
printf('  Number of examples: %d', L);

[dist,ratio] = dictdist(Dksvd,D);
printf('  Ratio of recovered atoms: %.2f%%\n', ratio*100);
