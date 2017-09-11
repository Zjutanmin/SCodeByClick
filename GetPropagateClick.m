% % % % 10000->1000
% % % % idx%%%%   10000*1
% % % % OriginalClickVector1   NQ*10000
% % % % 
% % % % Original_Index_IDX5000 = [];
% % % % ClusterRange = zeros(1, 1000);
% % % % for i = 1:1000
% % % %     tmp = find(idx== i);
% % % %     Original_Index_IDX5000 =[Original_Index_IDX5000; tmp];
% % % %      ClusterRange(i) = length(tmp);
% % % % end



function ClickPropagatedVector = GetPropagateClick( OriginalClickVector1, ...
    alpha, ClusterRange, Original_Index_IDX5000)
%FUN_VISUASIMILARITY�� ��������query��Ӧ���Ӿ��������֮����Ӿ����ƶ�
% Q1, Q2            :  �����query���ַ�����ʽ��
%Visual_Simi        :  �Ӿ����ƶ�

%(1)�ҵ����query��Ӧ���±�ţ����Q1��Q2��Query�±��ʾ������Ҫ��һ����

% if nargin < 1
%     OriginalClickVector1 = (rand(20, 484829));
% end
% if nargin < 2
%     alpha = 0.1;
% end
% if nargin < 3
%     load('Original_Index_IDX5000');
% end
% if nargin < 4
%     load('ClusterRange');
% end

% OriginalClickVector1 = sparse(rand(20, 484829));
%(2)����ʼ�ĵ����������cluster�ɿ����
% global Original_Index_IDX5000;
ClickPropagatedVector = OriginalClickVector1(:, Original_Index_IDX5000);

    
%(3)�������
% ClickPropagatedVector1 = Sub1_Fun_VisuaSimilarity(OriginalClickVector1, ClusterRange, alpha);
% 
% end
% 
% function ClickPropagatedVector = Sub1_Fun_VisuaSimilarity(xx, ClusterRange, alpha)
% %Sub1_Fun_VisuaSimilarity:����������������е������
% %xx                      :���ܵĴ������ĵ������
% %ClickPropagatedVector   ���ѵ������ĵ������
% 
% % global ClusterRange;
% ClickPropagatedVector = xx;

c_ClickPropagatedVector = cellfun(@(x) ClickPropagatedVector(:, x(1):x(2)), ...
	ClusterRange, 'UniformOutput', false);
L = cellfun(@(x) x(2)-x(1), ...
	ClusterRange, 'UniformOutput', false);
% c_ClickPropagatedVector = cellfun(@(x,y) (bsxfun(@plus, x, sum(x*alpha,2)/y)-x*alpha*(1/y+1)),...
% 	c_ClickPropagatedVector, L, 'UniformOutput', false);
for i = 1:length(L)
    if L{i}>0
        temp = c_ClickPropagatedVector{i} * alpha;
        c_ClickPropagatedVector{i} = bsxfun(@plus,c_ClickPropagatedVector{i},sum(temp,2)/L{i})-temp*(1/L{i}+1);
    end
end
ClickPropagatedVector = cell2mat(c_ClickPropagatedVector);
% for i = 1:length(ClusterRange)  %5000->10
%     r = ClusterRange{i};
%     startp = r(1);
%     endp = r(2);
%     block = xx(:, startp:endp);
%     
%     L = size(block, 2) - 1;
%     delta = sum(block*alpha,2)/L;
%     ClickPropagatedVector(:, startp:endp) = bsxfun(@plus, block, delta)-block*alpha*(1/L+1);
% end
end