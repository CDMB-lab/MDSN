function [mutualinfo] = mutualInfoTool(num,order,pts,class,values)
%MUTUALINFOTOOL �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
for j=1:num
    values(j,order+1)=MutualInformation(pts,class,values(j,1:order));
end
mutualinfo = values(:, order+1);
end

