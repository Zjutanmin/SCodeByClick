function Keachone = getKeachclass_img(type, K, Neachone)
%        Neachone ÿһ���ͼƬ��
%        type  ����ķ���
%        K ����� �� ÿ���ͼƬ��

    switch type
        case 'N'%%%ÿ�������һ��
            Keachone = min(ones(length(Neachone), 1)*K, Neachone);
        case 'P'%%%ÿ��ͼƬ��һ��
            Keachone = ceil(Neachone/K);
        case 'R'
            Keachone = ceil(Neachone * K);
    end
end