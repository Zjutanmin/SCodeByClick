function Keachone = getKeachclass_query(type, K, Neachone)
%        Neachone ÿһ���query��
%        type  ����ķ���
%        K ����� �� ÿ���query��

    switch type
        case 'N'%%%ÿ�������һ��
            Keachone = min(ones(length(Neachone), 1)*K, Neachone);
        case 'P'%%%ÿ��query��һ��
            Keachone = ceil(Neachone/K);
% %         case 'R'%%%ÿ��ķ������ɱ���
% %             Keachone = ceil(Neachone * K);
    end
end