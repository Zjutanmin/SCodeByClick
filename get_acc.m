function acc = get_acc(method,theta,alpha,P,i_Distype,cluster_method, type, T, gamma, textweight,normweight,SetUsePGNc,k_cluster_Q)
%      method         = method for determing cluster number in each category
%         N           : cluster number
%         P           : sample number in each cluster
%      theta          = parameter controlling cluster approach, when N>theta,sparse coding; while N<theta, K-means
%      alpha          = propagate ratio
%      P              = parameter for determing dictionary number n = gamma*N/P
%      i_Distype        = similarity functions on distance
%         1           : Sim = exp(-dist)
%         2           : Sim = 1./dist
%         3           : dist = mapminmax(dist(:),0,1)  & sim = exp(-dist)
%         4           : dist = mapminmax(dist(:),0,1)  & sim = 1./dist
%      cluster_method = query clustering method
%        1            : K-means
%        2            : hot query+sparse coding 
%        3            : KSVD+sparse coding
%      type           = propagation type
%        'WP'         : weighted propagation 
%        'EP'         : evenly propagation
%      T              = number of zero-item in sparse code coefficiences

%filename = ['IDX-',num2str(database_num),'-',num2str(cluster_num),'-',num2str(alpha),'-',num2str(net_num),'.mat'];
k_cluster_I = 20;
if nargin < 13
    k_cluster_Q = 20;
end
type_empty_cluster = 1;
global setting;
setting.type = i_Distype;

switch abs(cluster_method)
    case 1
        sname = 'KM';DefaultUsePGNc = 0;
    case 2
        sname = 'HS';DefaultUsePGNc = 1;
    case 3
        sname = 'KS';DefaultUsePGNc = 1;
end
% if SetUsePGNc
   
% else
%    
% end
if nargin < 12
    UsePGNc = DefaultUsePGNc;
else
    UsePGNc = SetUsePGNc;
end

setting.UsePGNc = UsePGNc;

if sign(cluster_method) == -1
    sname = [sname, method,'ME' ];
end
sname = [sname, method,'-' , type, '-k-', num2str(k_cluster_I), '-', num2str(k_cluster_Q)];
if alpha ~= 1
    sname = [sname, '-', num2str(alpha)];
end
if abs(cluster_method) ~= 1
    sname = [sname, 'ga-', num2str(gamma),...
        'T-', num2str(T),'th-', num2str(theta),'P-', num2str(P)];
end

sname = [sname,'D-', num2str(i_Distype)];
if abs(cluster_method)~=1
    sname = [sname, 'W-', num2str(textweight)];
end
if normweight~=1
    sname = [sname, 'N-', num2str(normweight)];
end

if UsePGNc~=DefaultUsePGNc
    sname = [sname, 'Nc-', num2str(UsePGNc)];
end

sname = [sname];


global isFortest;
try
    load(fullfile('result', [sname, 'c-', num2str(type_empty_cluster), '-acc.mat']), 'acc');       
%     acc
catch
    if isFortest == 1
        acc = 0;
    else
  

    IDX = mergeQ(sname, method,theta,alpha,P,i_Distype,cluster_method, type, T, gamma,textweight, normweight, k_cluster_Q);
    addpath('data')
    load('image_click_Dog283_0_click_non1_ND_fdatabase.mat', 'fdatabase1')%
    load('image_click_Dog283_0_RoundInfo3_1R_0.5.mat')%
    clear Valtset; %
    load('cluster_index.mat')%283
    load('query_index.mat')
    load('norm_Fea.mat')%

    image_index = [];

    for i = 1 : length(fdatabase1.label)   
        if(find(cluster_index == fdatabase1.label(i,1)))
            image_index(end+1,1) = i;            
        end
    end
    
    img_Fea_Clickcount(image_index,:) = 0;
    for i = 1 : size(img_Fea_Clickcount,2)
        if(nnz(img_Fea_Clickcount(:,i))==0)
        end
    end
    %save('query_index.mat','query_index');
    load('norm_Fea.mat')
    
    %total_label = fdatabase1.label(image_index);
    train_index = intersect(Trainset,image_index);
    train_label = fdatabase1.label(train_index);
 
    test_label = fdatabase1.label(intersect(Testset,image_index));



    %
    Train_Dataset = img_Fea_Clickcount(intersect(Trainset,image_index),query_index);
    Test_Dataset = img_Fea_Clickcount(intersect(Testset,image_index),query_index);
    %
    Train_Fea = [];Test_Fea = [];
    for i = 1 : max(IDX)
        index_class_num = find(IDX == i);
        Train_Fea = [Train_Fea,full(sum(Train_Dataset(:,index_class_num),2))];
        Test_Fea = [Test_Fea,full(sum(Test_Dataset(:,index_class_num),2))];
    end
    
    %
    Dis = EuDist2(Test_Fea,Train_Fea);

    
    [~,idx] = min(Dis,[],2);
    predict_label = train_label(idx);


    acc = 1 - nnz(predict_label - test_label)/ length(idx);  
    save(fullfile('result', [sname, 'c-', num2str(type_empty_cluster), '-acc.mat']), 'acc');
    end
end
end


