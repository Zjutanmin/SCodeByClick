
load('nonE_query_idx.mat');
load('Dog_ClickFea.mat','Dog_QQ_N');


load('N-C-k-10-20-cluster1-0.mat');
nonE_query_cluster = unique(IDX(nonE_query_idx));

s = 'nonE_textquery';
fid = fopen([s '.txt'],'wt');
t = [];

num_q = zeros(length(nonE_query_cluster),2);
for i = 1:length(nonE_query_cluster)
% %     fprintf(fid,'The cluster of query: %s\n',num2str(nonE_query_cluster(i)));
    q = find(IDX==nonE_query_cluster(i));
    nonE_query = nonE_query_idx(find(IDX(nonE_query_idx)==nonE_query_cluster(i)));
    E_query = setdiff(q,nonE_query);
    if isempty(E_query)
        t = [t;nonE_query];
        num_q(i,:) = [length(nonE_query),length(E_query)];
        r = randperm(length(E_query),min(length(E_query),1));
% % %         fprintf(fid,'non E (#%s) / ',num2str(length(nonE_query)));
% % %         fprintf(fid,'E (#%s):+ ',num2str(length(E_query)));
% % %         fprintf(fid,'%s/%s\n\n',num2str(length(nonE_query)),num2str(length(q)));
        fprintf(fid,'(%d,%d,%.2f%%)\n',length(nonE_query),length(E_query),length(nonE_query)/length(q)*100);
        for i_nonE_query = 1:length(nonE_query)
            fprintf(fid,'%s\n',Dog_QQ_N{nonE_query(i_nonE_query)});
        end
        if ~isempty(r)
% % %             fprintf(fid,'..................................................\n');
% % %             for i_E_query = 1:length(r)
% % %                 fprintf(fid,'%s\n',Dog_QQ_N{E_query(i_E_query)});
% % %             end
            fprintf(fid,'%s, ...\n',Dog_QQ_N{E_query(r)});
        else
            fprintf(fid,'\n');
        end
        fprintf(fid,'--------------------------------------------------\n');
    end
end
save(s,'num_q');
clear IDX;
fclose(fid);

nonE_query_idx = t;



load('result/N-C-k-10-20-cluster1.mat');
nonE_query_cluster = unique(IDX(nonE_query_idx));

s = 'nonE_clickquery';
fid = fopen([s '.txt'],'wt');

num_q = zeros(length(nonE_query_cluster),2);
for i = 1:length(nonE_query_cluster)
% %     fprintf(fid,'The cluster of query: %s\n',num2str(nonE_query_cluster(i)));
    q = find(IDX==nonE_query_cluster(i));
    nonE_query = nonE_query_idx(find(IDX(nonE_query_idx)==nonE_query_cluster(i)));
    E_query = setdiff(q,nonE_query);
    num_q(i,:) = [length(nonE_query),length(E_query)];
    r = randperm(length(E_query),min(length(E_query),1));
% % %     fprintf(fid,'non E (#%s) / ',num2str(length(nonE_query)));
% % %     fprintf(fid,'E (#%s): ',num2str(length(E_query)));
% % %     fprintf(fid,'%s/%s\n\n',num2str(length(nonE_query)),num2str(length(q)));
    fprintf(fid,'(%d,%d,%.2f%%)\n',length(nonE_query),length(E_query),length(nonE_query)/length(q)*100);
    for i_nonE_query = 1:length(nonE_query)
        fprintf(fid,'%s\n',Dog_QQ_N{nonE_query(i_nonE_query)});
    end
    if ~isempty(r)
% % %             fprintf(fid,'..................................................\n');
% % %             for i_E_query = 1:length(r)
% % %                 fprintf(fid,'%s\n',Dog_QQ_N{E_query(i_E_query)});
% % %             end
        fprintf(fid,'%s, ...\n',Dog_QQ_N{E_query(r)});
    else
        fprintf(fid,'\n');
    end
    fprintf(fid,'--------------------------------------------------\n');

end
save(s,'num_q');
clear IDX;
fclose(fid);



load('N-C-k-10-20-cluster1-0.001.mat');
nonE_query_cluster = unique(IDX(nonE_query_idx));

s = 'nonE_mixquery';
fid = fopen([s '.txt'],'wt');

num_q = zeros(length(nonE_query_cluster),2);
for i = 1:length(nonE_query_cluster)
% %     fprintf(fid,'The cluster of query: %s\n',num2str(nonE_query_cluster(i)));
    q = find(IDX==nonE_query_cluster(i));
    nonE_query = nonE_query_idx(find(IDX(nonE_query_idx)==nonE_query_cluster(i)));
    E_query = setdiff(q,nonE_query);
    num_q(i,:) = [length(nonE_query),length(E_query)];
    r = randperm(length(E_query),min(length(E_query),1));
% % %     fprintf(fid,'non E (#%s) / ',num2str(length(nonE_query)));
% % %     fprintf(fid,'E (#%s): ',num2str(length(E_query)));
% % %     fprintf(fid,'%s/%s\n\n',num2str(length(nonE_query)),num2str(length(q)));
    fprintf(fid,'(%d,%d,%.2f%%)\n',length(nonE_query),length(E_query),length(nonE_query)/length(q)*100);
    for i_nonE_query = 1:length(nonE_query)
        fprintf(fid,'%s\n',Dog_QQ_N{nonE_query(i_nonE_query)});
    end
    if ~isempty(r)
% % %             fprintf(fid,'..................................................\n');
% % %             for i_E_query = 1:length(r)
% % %                 fprintf(fid,'%s\n',Dog_QQ_N{E_query(i_E_query)});
% % %             end
        fprintf(fid,'%s, ...\n',Dog_QQ_N{E_query(r)});
    else
        fprintf(fid,'\n');
    end
    fprintf(fid,'--------------------------------------------------\n');

end
save(s,'num_q');
clear IDX;
fclose(fid);



