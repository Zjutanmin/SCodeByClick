function Keachone = getKeachclass_img(type, K, Neachone)
%        Neachone 每一类的图片数
%        type  聚类的方法
%        K 类别数 或 每类的图片数

    switch type
        case 'N'%%%每类类别数一样
            Keachone = min(ones(length(Neachone), 1)*K, Neachone);
        case 'P'%%%每类图片数一样
            Keachone = ceil(Neachone/K);
        case 'R'
            Keachone = ceil(Neachone * K);
    end
end