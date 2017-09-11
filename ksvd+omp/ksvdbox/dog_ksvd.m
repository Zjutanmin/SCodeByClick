addpath('D:\m\code\ksvd+omp\ksvdbox')
addpath('D:\m\code\ksvd+omp\ompbox')
load('Dict.mat');
load('image_click_Dog283_0_img_Fea_Clickcount.mat')

% % X = img_Fea_Clickcount';

X = img_Fea_Clickcount';

k = 4;
m = size(D,1);
params.lamdaRange = [0.1:0.1:1];

params.data = X;           % %���ݱ���
params.Tdata = k;           % %Լ��������Լ��ֵ
params.dictsize = m;        % %�ֵ���Ҫ�������ĸ���(e1,e2.....em)
params.iternum = 30;        % % ��������
params.memusage = 'high';
params.fixedcode = 0;          % %���벻��
params.method_dic = 'L0';      % %����0���Ƿ���1�����ַ�ʽ
params.SCode.Code = [];  %%%Min Tan
params.lamda = 1;         % %��Ҫ���ڣ�����Ϊ1ʱ��Ĳ���

params.fixedDic =0;
[~,code] = ksvd(params,'');    %%%%���

%% show results %%

figure; plot(err); title('K-SVD error convergence');
xlabel('Iteration'); ylabel('RMSE');

printf('  Dictionary size: %d x %d', n, m);
printf('  Number of examples: %d', L);

[dist,ratio] = dictdist(Dksvd,D);
printf('  Ratio of recovered atoms: %.2f%%\n', ratio*100);
