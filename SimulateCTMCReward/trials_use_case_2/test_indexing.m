clear; clc;
a = rand(1e5,1);
b = rand(1e5,1);
index_1 = a>.5;
addpath('..\CopyMask')
%% Test
tic;
c = zeros(1e5,1);
for i = 1:1e4
    c = a;
end
toc

tic;
d = zeros(sum(index_1),1);
for i = 1:1e4
    d = a(index_1);
end
toc

tic;
for i = 1:1e4
    b = a+b;
end
toc

tic;
for i = 1:1e4
    b(index_1) = a(index_1)+b(index_1);
end
toc

tic;
for i = 1:1e4
    b(index_1) = CopyMask(a,index_1)+CopyMask(b,index_1);
end
toc

tic;
for i = 1:1e4
    b = a;
end
toc

tic;
for i = 1:1e4
    b(index_1) = a(index_1);
end
toc

tic;
for i = 1:1e4
    b(1:5000) = b(1:5000).*a(1:5000);
end
toc

tic
for i = 1:1e4
    b = b.*a;
end
toc
