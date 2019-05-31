function [x,y,z] = find3(X)
x = [];
y = [];
z = [];
% a = size(X,1)*size(X,2);
for i = 1:size (X,3)
    temp = X(:,:,i);
    [x_temp,y_temp] = find(temp);
    z_temp = i*ones(size(x_temp));
    x = [x;x_temp];
    y = [y;y_temp];
    z = [z;z_temp];
end