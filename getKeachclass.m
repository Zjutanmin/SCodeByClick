function Keachone = getKeachclass(type, K, Neachone)
%        Neachone ÿһ���ͼƬ��
%        type  ����ķ���
%        K ����� �� ÿ���ͼƬ��


    switch type
        case 'equalN'%%%ÿ�������һ��
            Keachone = min(ones(length(Neachone), 1)*K, Neachone);
        case 'equalPN'%%%ÿ��ͼƬ��һ��
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