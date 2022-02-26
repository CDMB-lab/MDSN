function [mutualinfo] = mutualInfoTool(num,order,pts,class,values)
%MUTUALINFOTOOL 此处显示有关此函数的摘要
%   此处显示详细说明
for j=1:num
    values(j,order+1)=MutualInformation(pts,class,values(j,1:order));
end
mutualinfo = values(:, order+1);
end

