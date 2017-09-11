global isFortest;
if ~exist('istest', 'var')
    isFortest = 0;
else
    isFortest = istest;
end
Rresult = [];
% % % % %%%%\gamma where
% % % % %1.propagated click feature with varying \alpha
% % % % mkdir('result');
theta = Inf;method = 'N';P = 1;cluster_method = 1;gamma = 0; 
i_alpha = 0;i_Distype = 1;
i = 0;acc = [];
i = i+1;acc(i) = get_acc(method,theta,i_alpha,P,i_Distype,cluster_method, 'O',1,gamma,1,2,0,10);
i = i+1;acc(i) = get_acc(method,theta,i_alpha,P,i_Distype,cluster_method, 'C',1,gamma,1,2,0,10);
Rresult = acc;



%%%
theta = Inf;method = 'N';P = 1;cluster_method = 1;gamma = 0; 
alpha = [-0.9:0.2:0.9];

Distype = [1];

i = 0;acc = [];
for i_Distype = Distype
for i_alpha = alpha
%     i_alpha
    i = i+1;acc(i) = get_acc(method,theta,i_alpha,P,i_Distype,cluster_method, 'P',1,gamma,1,2,0,10);
    i = i+1;acc(i) = get_acc(method,theta,i_alpha,P,i_Distype,cluster_method, 'WP',1,gamma,1,2,0,10);
end
end
res = reshape(acc, 2, []);
res = reshape(res, [2, 5, 2]);

RES1 = res(1,:, 1)*100; %%%P -
RES1 = RES1(end:-1:1);
RES2 = res(1,:, 2)*100; %%%p +
RES3 = res(2,:, 1)*100; %%%wP -
RES3 = RES3(end:-1:1);
RES4 = res(2,:, 2)*100; %%%wp +
plot(alpha(6:end), RES1, 'r-','linewidth', 2)



hold on;plot(alpha(6:end), RES2,  'go-','linewidth', 2)
hold on;plot(alpha(6:end), RES3, 'b-','linewidth', 2)
hold on;plot(alpha(6:end), RES4, 'k-','linewidth', 2)
legend({'P-', 'P+', 'WP-', 'WP+'}, 'fontsize', 16);
xlabel('\alpha' ,'fontsize', 16)
ylabel('Accuracy (%)', 'fontsize', 16)
result = reshape(acc(1:2:end), [length(alpha), length(Distype)]);
[Mresult, ind] = max(result(:));[x,y] = ind2sub(size(result), ind);
Rresult = [Rresult, Mresult];


result = reshape(acc(2:2:end), [length(alpha), length(Distype)]);
[Mresult, ind] = max(result(:));[x,y] = ind2sub(size(result), ind);
Rresult = [Rresult, Mresult];


alpha = alpha(x);Distype = Distype(y);


Rresult1 = Rresult;

%%%上面已经跑了





%%2.P and T for hot-query + sparse coding                跑完这个选P，T
theta = 600;method = 'N';cluster_method = [-2];normweight = [2]; Distype = [1];           
P = [1, 5,10,15];
T = [0.1:0.05:0.3];textweight = [1];gamma =[0.1]; alpha = -0.5;
i = 0;acc = [];
for i_P = P
%     i_P
    for i_T = T
%         i_T        
            for i_textweight = textweight%%%%FOR BIGGER
                for i_normweight = normweight
                i = i+1;
                acc(i) = get_acc(method,theta,alpha,i_P,Distype,-3, 'WP',i_T, gamma, i_textweight,i_normweight,0, 10);
                i = i+1;
                acc(i) = get_acc(method,theta,alpha,i_P,Distype,-2, 'WP',i_T, gamma, i_textweight,i_normweight,0, 10);
%                i = i+1;
%                acc(i) = get_acc(method,theta,alpha,i_P,Distype,2, 'WP',i_T, gamma, i_textweight,i_normweight,0, 10);
                end           
        end
    end
end
res = reshape(acc,2,[])

Nlen = length(T)* length(gamma)* length(textweight)* length(normweight)*2;

%%%%For P = 1
result = reshape(acc(1:2:Nlen), [length(normweight),length(textweight), length(gamma), length(T), 1]);
[Mresult, ind] = max(result(:));[u, x,y,z,w] = ind2sub(size(result), ind);
Rresult = [Rresult, Mresult];


P1 = P(1);
normweight_SVD_S = normweight(u);textweight_SVD_S = textweight(x);gamma_SVD_S = gamma(y);T_SVD_S = T(z);P_SVD_S = P1(w);

result = reshape(acc(2:2:Nlen), [length(normweight),length(textweight), length(gamma), length(T), 1]);
[Mresult, ind] = max(result(:));[u, x,y,z,w] = ind2sub(size(result), ind);
Rresult = [Rresult, Mresult];



P1 = P(1);
normweight_HOT_S = normweight(u);textweight_HOT_S = textweight(x);gamma_HOT_S = gamma(y);T_HOT_S = T(z);P_HOT_S = P1(w);



%%% For P > 1
result = reshape(acc(Nlen+1:2:end), [length(normweight),length(textweight), length(gamma), length(T), length(P)-1]);
[Mresult, ind] = max(result(:));[u, x,y,z,w] = ind2sub(size(result), ind);
Rresult = [Rresult, Mresult];



P1 = P(2:end);
normweight_SVD_M = normweight(u);textweight_SVD_M = textweight(x);gamma_SVD_M = gamma(y);T_SVD_M = T(z);P_SVD_M = P1(w);


result = reshape(acc(Nlen+2:2:end), [length(normweight),length(textweight), length(gamma), length(T), length(P)-1]);
[Mresult, ind] = max(result(:));[u, x,y,z,w] = ind2sub(size(result), ind);
Rresult = [Rresult, Mresult];


P1 = P(2:end);
normweight_HOT_M = normweight(u);textweight_HOT_M = textweight(x);gamma_HOT_M = gamma(y);T_HOT_M = T(z);P_HOT_M = P1(w);

T_SVD_S = 0.2;
T_HOT_S = 0.2;
T_SVD_M = 0.2;
T_HOT_M = 0.2;
P_SVD_S = 1;
P_HOT_S = 1;
P_SVD_M = 10;
P_HOT_M = 10;

%%%给我看完结果后，开始跑下面的
%4.overall performance            10 
alpha = -0.5;          
normweight = 2;textweight  = 1;gamma = 0.1;
method = 'N';Distype = [1];type = 'WP';
theta = [100:100:1000];
i = 0;acc = [];
for i_theta = theta
    i_theta
    i = i+1;acc(i) = get_acc(method,i_theta,alpha,P_SVD_S,Distype,-3, 'WP',T_SVD_S,gamma,textweight,1,0, 10);
    i = i+1;acc(i) = get_acc(method,i_theta,alpha,P_HOT_S,Distype,-2, 'WP',T_HOT_S,gamma,textweight,1,0, 10);
end

res = reshape(acc,2,[])
res1 = res;
result = acc(1:2:end);
[Mresult, ~] = max(result(:));
Rresult = [Rresult, Mresult];

result = acc(2:2:end);
[Mresult, ~] = max(result(:));
Rresult = [Rresult, Mresult];



%%%%FOR part 3, PLEASE SET ,alpha, Distype ,textweight, gamma ,T and P

%3.KSVD / HOT under different \theta             10 
% Rresult = [];
alpha = -0.5;          
normweight = 2;textweight  = 1;gamma = 0.1;
method = 'N';Distype = [1];type = 'WP';
theta = [100:100:1000];
i = 0;acc = [];
for i_theta = theta
    i_theta
    i = i+1;acc(i) = get_acc(method,i_theta,alpha,P_SVD_M,Distype,-3, 'WP',T_SVD_M,gamma,textweight,1,0, 10);
    i = i+1;acc(i) = get_acc(method,i_theta,alpha,P_HOT_M,Distype,-2, 'WP',T_HOT_M,gamma,textweight,1,0, 10);
end

res = reshape(acc,2,[])
res2 = res;
result = acc(1:2:end);
[Mresult, ~] = max(result(:));
Rresult = [Rresult, Mresult];

result = acc(2:2:end);
[Mresult, ~] = max(result(:));
Rresult = [Rresult, Mresult];