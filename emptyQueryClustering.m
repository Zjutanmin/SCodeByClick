function [C, A] = emptyQueryClustering(QueryFea, query_g1_empty_idx,Keachone_query,len,gamma,P,normweight)
    Keachone = round(mean(Keachone_query));
    data = (QueryFea');
    number = min(Keachone*len,length(query_g1_empty_idx));
    data = Mynomalize(data, 1, normweight);
    [C,A] = vl_kmeans(full(data),number);
%     [C,A] =  Myclustering(data,number,1,1,gamma,P,0,0.1,normweight);
%            % Myclustering(temp_QueryFea,Keachone_query,i,cluster_method,gamma,P,theta,T, textweight,normweight)   
    
    C = C';
    A = A';
end