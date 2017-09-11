method_A = {'N', 'P'};
type_A = {'C', 'P'};
c_k_cluster_I_T = {[10, 20, 50, 100], [5, 10, 50, 100, 500, 1000]};
c_k_cluster_Q_T = {[10, 20, 50, 100], [5, 10, 20, 50, 500, 1000]}; 

mkdir('query_cluster_text');

load('Dog_ClickFea.mat');

sname = [method_A{1} '-' type_A{1} '-k-'];
k_cluster_I = c_k_cluster_I_T{1};
k_cluster_Q = c_k_cluster_Q_T{1};
path = pwd;


for i_k_cluster_I = 1:length(k_cluster_I)
    for i_k_cluster_Q = 1:length(k_cluster_Q)
        s = [sname num2str(k_cluster_I(i_k_cluster_I)) '-' num2str(k_cluster_Q(i_k_cluster_Q)) '-cluster1'];
        load(['click_idx/' s '.mat'],'IDX');
        s_1 = fullfile('query_cluster_text', s);
        mkdir(s_1);
        cd(s_1);
        for i = 1:max(IDX)
            fid = fopen([num2str(i) '.txt'],'wt');
            fprintf(fid,'%s\n',Dog_QQ_N{IDX==i});
            fclose(fid);
        end
        cd(path);
        disp(s);
    end
end
