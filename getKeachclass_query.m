function Keachone = getKeachclass_query(type, K, Neachone)
%        Neachone 每一类的query数
%        type  聚类的方法
%        K 类别数 或 每类的query数

    switch type
        case 'N'%%%每类类别数一样
            Keachone = min(ones(length(Neachone), 1)*K, Neachone);
        case 'P'%%%每类query数一样
            Keachone = ceil(Neachone/K);
% %         case 'R'%%%每类的分类数成比例
% %             Keachone = ceil(Neachone * K);
    end
end