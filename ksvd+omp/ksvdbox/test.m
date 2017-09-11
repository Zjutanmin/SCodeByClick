data = rand(112,231);tempcount = sum(data,1);
pX = data';pk = 15;
params.data = pX;
pk = 6;
       
Nc = 6;hotrate = 0.15;
[IDX,C ] = kmeans(params.data',Nc); 
            [a,b,c] = unique(IDX);
            pD = []; code_size = [];        
            
             for pi = 1:length(a)
                 index = find(c == pi);
                 [~, iind] = sort(tempcount( index), 'descend');
                 iind = iind(1:max(1, floor(length(index)*hotrate)));
                 pD = [pD,params.data(:, index(iind)) ]; 
                 code_size(pi) = length(iind);
                 code_size1(pi) = size(pD,2);
             end
             code_size2(1) = 0;
             for pi = 2:length(code_size1)+1
                 code_size2(pi) = code_size1(pi-1);
             end
            
             params.Tdata = Nc;
             params.dictsize = size(pD,2);
             params.initdict = pD;
                 [~ ,code] = ksvd(params,'');
                code1 = [];
                 
              for pi = 1:length(a)     
                  start1 = 1+code_size2(pi);
                  end1 = code_size2(pi+1);
                  code1 = [code1;sum(      code(start1:end1,:)  ,1)];
                  
              end
              
             [~,cluster_idx] = max(code1, [], 1);     