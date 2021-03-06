% clear,clc
filename = 'bike2.txt';
delimiterIn = ' ';
headerlinesIn = 1;
A = importdata(filename,delimiterIn,headerlinesIn);

X=A.data;

obst_i=[];
for i=1:0.5:4
    obst_i=[obst_i; [i,4]];
end
for i=1:0.5:4
    obst_i=[obst_i; [1,i]];
end
for i=-1:-0.5:-4
    obst_i=[obst_i; [i,4]];
end
for i=1:0.5:4
    obst_i=[obst_i; [-1,i]];
end
        
figure(1)
hold on
plot(X(:,1),X(:,2))
plot(obst_i(:,1), obst_i(:,2), '*r')
grid on


figure(2)
hold on
plot(X(:,4))
plot(X(:,5))
legend('V','gam')
grid on