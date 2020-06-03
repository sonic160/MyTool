clear; clc;
A = rand(5,1e5);
for i = 1:1e5
    A(1,:) = rand(1,1e5);
end
B = rand(1e5,5);
for i = 1:1e5
    B(:,1) = rand(1e5,1);
end