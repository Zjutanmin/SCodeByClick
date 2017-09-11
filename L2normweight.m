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