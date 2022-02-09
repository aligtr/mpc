% clear,clc
filename = 'cppstudio.txt';
delimiterIn = ' ';
headerlinesIn = 1;
A = importdata(filename,delimiterIn,headerlinesIn);

X=A.data;
hold on
plot(X(:,1),X(:,2))