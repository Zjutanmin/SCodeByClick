function ksvddemo
%KSVDDEMO K-SVD training demonstration.
%  KSVDDEMO generates random sparse examples over a random overcomplete
%  dictionary, adds noise, and uses K-SVD to recover the dictionary from
%  the examples.
%
%  To run the demo, type KSVDDEMO from the Matlab prompt.
%
%  See also KSVDDENOISEDEMO.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  May 2009


disp(' ');
disp('  **********  K-SVD Demo  **********');
disp(' ');
disp('  This demo generates a random dictionary and random sparse examples over this');
disp('  dictionary. It then adds noise to the examples, and uses K-SVD to recover');
disp('  the dictionary. The demo plots the convergence of the K-SVD target function,');
disp('  and computes the fraction of correctly recovered atoms.');
disp(' ');
load('image_click_Dog283_0_img_Fea_Clickcount.mat')

% dictionary dimensions  ά��
n = 20;
m = 50;

% number of examples
L = 1500;

% sparsity of each example   ��3ά��ʾ
k = 3;
% noise power (dB)
snr = 20;



%% generate random dictionary and data %%   �ֵ�����

D = normcols(randn(n,m));

Gamma = zeros(m,L);
for i = 1:L
  p = randperm(m);
  Gamma(p(1:k),i) = randn(k,1);
end

X = D*Gamma;

X = normcols(X) + 10^(-snr/20)*normcols(randn(n,L));




%% run k-svd training %%
% % load('Dict.mat')
% % load('image_click_Dog283_0_img_Fea_Clickcount.mat')
% % X = img_Fea_Clickcount;
% % addpath('D:\m\code\ksvd+omp\ksvdbox')
% % addpath('D:\m\code\ksvd+omp\ompbox')
params.lamdaRange = [0.1:0.1:1];
addpath('D:\m\code\ksvd+omp\ompbox');
params.data = X;           % %���ݱ���
params.Tdata = k;           % %Լ��������Լ��ֵ
params.dictsize = m;        % %�ֵ���Ҫ�������ĸ���(e1,e2.....em)
params.iternum = 30;        % % ��������
params.memusage = 'high';
params.fixedcode = 0;          % %���벻��
params.method_dic = 'L0';      % %����0���Ƿ���1�����ַ�ʽ
params.SCode.Code = [];  %%%Min Tan
params.lamda = 1;         % %��Ҫ���ڣ�����Ϊ1ʱ��Ĳ���
[Dksvd,g,err] = ksvd(params,'');    %%%%���


%% show results %%

figure; plot(err); title('K-SVD error convergence');
xlabel('Iteration'); ylabel('RMSE');

printf('  Dictionary size: %d x %d', n, m);
printf('  Number of examples: %d', L);

[dist,ratio] = dictdist(Dksvd,D);
printf('  Ratio of recovered atoms: %.2f%%\n', ratio*100);
