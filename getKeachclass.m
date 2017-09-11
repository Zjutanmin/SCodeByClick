function Keachone = getKeachclass(type, K, Neachone)
%        Neachone 每一类的图片数
%        type  聚类的方法
%        K 类别数 或 每类的图片数


    switch type
        case 'equalN'%%%每类类别数一样
            Keachone = min(ones(length(Neachone), 1)*K, Neachone);
        case 'equalPN'%%%每类图片数一样
            Keachone = round(Neachone/K);
    end
    
    load('image_click_Dog283_0_CNN_Alex1_ND_S_S1_data.mat');
    
%     for i = 1:max(img_label) 
%         index = find(data_label == i);
%         fea = data_fea(index,:)';
%         [C, A] = vl_kmeans(fea,Keachone(i));
%     end
    
    [C,A] = arrayfun(@(x) vl_kmeans(data_fea(find(data_label==x),:)',Keachone(x)),1:max(data_label),'UniformOutput',false);
    
    
    
end