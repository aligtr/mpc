% clear,clc
filename = 'dyn1.txt';
delimiterIn = ' ';
headerlinesIn = 1;
A = importdata(filename,delimiterIn,headerlinesIn);

X=A.data;

obst_i=[];
for i=1:0.5:5
    obst_i=[obst_i; [i,4]];
end
for i=1:0.5:4
    obst_i=[obst_i; [1,i]];
end
for i=-1:-0.5:-5
    obst_i=[obst_i; [i,4]];
end
for i=1:0.5:4
    obst_i=[obst_i; [-1,i]];
end

subplot(1,3,1)
hold on
plot(X(:,1),X(:,2))
plot(obst_i(:,1), obst_i(:,2), '*r')
grid on

subplot(1,3,2)
hold on
plot(X(:,3))
plot(X(:,4))
plot(X(:,5))
plot(X(:,6))
legend('M1','M2','M3','M4')
grid on

subplot(1,3,3)
hold on
plot(X(:,7))
plot(X(:,8))
plot(X(:,9))
plot(X(:,10))
legend('UG1','UG2','UG3','UG4')
grid on