theta = Inf;method = 'N';P = 1;cluster_method = 1;gamma = 0; 
alpha = [-0.9:0.2:0.9];alpha = alpha(1);
%alpha = [0.1:0.2:0.9];
Distype = [1];

i = 0;acc = [];
for i_Distype = Distype
for i_alpha = alpha
    i_alpha
    i = i+1;acc(i) = get_acc(method,theta,i_alpha,P,i_Distype,cluster_method, 'P',1,gamma,1,1,0,10);
    %i = i+1;acc(i) = get_acc(method,theta,i_alpha,P,i_Distype,cluster_method, 'WP',1,gamma,1,1,0,10);
end
end
res = reshape(acc, 2, []);