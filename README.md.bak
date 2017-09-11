# Query-Modeling-for-Click-Data





This is the code for the Click data based Query Merging alogrithm.



AUTHORS:
Min Tan <tanmin@hdu.edu.cn>


INTRODUCTION
------------
This package contains the MATLAB code for click feature exaction and hot-query based clustering.

The software included here implements the algorithm described in

[1] Weichen WU, Min Tan, Jun Yu, Guangjian Zheng. 
   Query Modeling for Click Data based Image Recognition using Graph based Propagation and Sparse Coding. 
   International Conference on Internet Multimedia Computing and Service 2017.




INSTALLATION
------------

1. Requirements

This software requires MATLAB R2007a or later.  Because it makes extensive use
of the "bsxfun" function, earlier versions of Matlab will not work.


2. Running the demo
Go to the demo files '/path' to run the demo script:

    >> demo




NOTE:
1).Related figure and results are saved in '/path\\result'
2).Related data are saved in '/path\\data'
3).The text feature for queries are saved in "/path/TextFea\image_click_Dog283_0_QueryFea"
4).The image clustering result is saved in "/path/ImgC"



3. Important functions
MAIN FUNCTION - COMPUTE ACCURACY
--------
acc = get_acc(method,theta,alpha,P,i_Distype,cluster_method, type, T, gamma, textweight,normweight,SetUsePGNc,k_cluster_Q)
%      method         = method for determing cluster number in each category
%         N           : cluster number
%         P           : sample number in each cluster
%      theta          = parameter controlling cluster approach, when N>theta,sparse coding; while N<theta, K-means
%      alpha          = propagate ratio
%      P              = parameter for determing dictionary size = N_k*P
%      i_Distype        = similarity functions on distance
%         1           : Sim = exp(-dist), default
%         2           : Sim = 1./dist
%         3           : dist = mapminmax(dist(:),0,1)  & sim = exp(-dist)
%         4           : dist = mapminmax(dist(:),0,1)  & sim = 1./dist
%      cluster_method = query clustering method
%        1            : K-means
%        2            : hot query + sparse coding 
%        3            : KSVD + sparse coding
%      type           = propagation type
%        'WP'         : weighted propagation 
%        'P'          : evenly propagation
%      T              = number of zero-item in sparse code coefficiences
%      gamma          = parameter for determing query cluster number
%      textweight     = weight for text feature respect to click feature in query representation for roughly grouping     
%      normweight     = normlization type
          1           : L1 normalization, default
%         2           : L2 normalization                
%      SetUsePGNc     = query cluster number assignment
          0           : N_k = k_cluster_Q
%         1                      : N_k = N*gamma/P
       k_cluster_Q   = pre-defined query cluster number, default value is 20
       
       acc            = Per-class accuracy



QUERY MERGING
--------
IDX = mergeQ(sname, method,theta,alpha,P,Distype,cluster_method,type, T,gamma, textweight,normweight,k_cluster_Q)
%      sname         = filename for storing results
      method         = method for determing cluster number in each category
%         N           : cluster number
%         P           : sample number in each cluster
%      theta          = parameter controlling cluster approach, when N>theta,sparse coding; while N<theta, K-means
%      alpha          = propagate ratio
%      P              = parameter for determing dictionary size = N_k*P
%      Distype        = similarity functions on distance
%         1           : Sim = exp(-dist), default
%         2           : Sim = 1./dist
%         3           : dist = mapminmax(dist(:),0,1)  & sim = exp(-dist)
%         4           : dist = mapminmax(dist(:),0,1)  & sim = 1./dist
%      cluster_method = query clustering method
%        1            : K-means
%        2            : hot query + sparse coding 
%        3            : KSVD + sparse coding
%      type           = propagation type
%        'WP'         : weighted propagation 
%        'P'          : evenly propagation
%      T              = number of zero-item in sparse code coefficiences
%      gamma          = parameter for determing query cluster number
%      textweight     = weight for text feature respect to click feature in query representation for roughly grouping     
%      normweight     = normlization type
          1           : L1 normalization, default
%         2           : L2 normalization                
       k_cluster_Q   = pre-defined query cluster number, default value is 20
       
       IDX            = cluster index for each query


CLICK PROPAGATION
--------
PropagatedVector_W = GetPropagate_W( OriginalVector, alpha, dist, clusterRange, Distype)
%      OriginalVector = Original click feature
%      alpha          = propagate ratio
%      dist           = distance matrix
%      Distype        = similarity functions on distance matrix
%         1           : Sim = exp(-dist), default
%         2           : Sim = 1./dist
%         3           : dist = mapminmax(dist(:),0,1)  & sim = exp(-dist)
%         4           : dist = mapminmax(dist(:),0,1)  & sim = 1./dist
       clusterRange   = cluster 
       
       PropagatedVector_W= propagated click feature

QUERY CLUSTERING BASED ON CLICK FEATURE
--------
[cluster_vector,cluster_idx_q] = Myclustering(temp_QueryFea,...
    Keachone_query,i,cluster_method,gamma,P,theta,T, textweight,normweight,ClickCount)   
%      temp_QueryFea  = propagated click feature for queries in one category
       Keachone_query = pre-defined query cluster number, default value is 20
%      cluster_method = query clustering method
%        1            : K-means
%        2            : hot query + sparse coding 
%        3            : KSVD + sparse coding
%      gamma          = parameter for determing query cluster number
%      P              = parameter for determing dictionary size = N_k*P
       theta          = parameter controlling cluster approach, when N>theta,sparse coding; while N<theta, K-means
%      T              = number of zero-item in sparse code coefficiences
%      textweight     = weight for text feature respect to click feature in query representation for roughly grouping     
%      normweight     = normlization type
          1           : L1 normalization, default
%         2           : L2 normalization                
       ClickCount     = total click count for queries in one category
       
       cluster_vector = cluster center vector
       cluster_idx_q  = cluster index for each queryPropagatedVector_W