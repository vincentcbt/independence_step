%dual extra 单次循环版
clear;clc;
tic
nodes = 100;
L = readmatrix('E:\My_Research\Programe\Distributed_optimization\graph_L.csv');
L(1,:)=[];
L_eig = eig(L);
lambda_m = max(L_eig);
W = eye(nodes) - L/lambda_m;
% W_tao = (eye(nodes) + W)/2;
load('D:\MATLAB_2022a\cbt_works\optx_100nods.mat')
load('D:\MATLAB_2022a\cbt_works\optx_100_cons.mat')
%随机生成参数
rng(123); %随机数种子
% a = (0.8-0.3).*rand(nodes,1)+0.3;  %生成[0.3,0.8]均匀分布随机数
% b = (10-1).*rand(nodes,1)+1;
% x_lb = (5-(-5)).*rand(nodes,1)-5;
% x_ub = x_lb + 5;

a = fix((0.3 + (0.8-0.3)*rand(nodes,1))*100)/100;
b = fix((1 + (10-1)*rand(nodes,1))*100)/100;
x_lb = fix(((5-(-5)).*rand(nodes,1)-5)*100)/100;
x_ub = x_lb+5;

% beta = 0.1;
beta = a;
% beta = 0.4;
cc = 1;
W_tao = eye(nodes) - cc*beta.*eye(nodes) + cc*beta.*W;

x0 = zeros(nodes,1);
mu0 = zeros(nodes,1);

mu_new1 = zeros(nodes,1);
mu_new2 = zeros(nodes,1);
x_1 = zeros(nodes,1);
x_new = zeros(nodes,1);
y0 = zeros(nodes,1);
y_new = zeros(nodes,1);

x_value = [];
mu_value = [];
err_value = [];

r = rand(1,100);
r = (100*r/sum(r))'*1; 
% r = ones(nodes,1);

%% dual extra 
parfor i=1:nodes
    fun = @(x) a(i)*x^2 - b(i)*x + mu0(i)*(x-r(i));
%     有约束
%             option = optimoptions("fmincon","Display","off");
%             opt_x = fmincon(fun,x0(i),[],[],[],[],x_lb(i),x_ub(i),[],option);
%             x_new(i) = opt_x;
%     无约束
    option = optimoptions("fmincon","Display","off");
    opt_x = fmincon(fun,x0(i),[],[],[],[],-3,3,[],option);
    x_1(i) = opt_x;
end


mu_new1 = W*mu0 + beta.*(x_1 - r);

k = 1;
while true
    parfor i=1:nodes
        fun = @(x) a(i)*x^2 - b(i)*x + mu_new1(i)*(x-r(i));
        option = optimoptions("fmincon","Display","off");
        opt_x = fmincon(fun,0,[],[],[],[],-3,3,[],option);
        x_new(i) = opt_x;
    end
    mu_new2 = W_tao*(2*mu_new1 - mu0 + beta.*(x_new - x_1));

    if norm(mu_new2-mu0,2) < 1e-6
        break
    else
        k = k+ 1;
        mu0 = mu_new1;
        mu_new1 = mu_new2;
        x_1 = x_new;
        x0 = x_1;
%         err = norm(x0-x_star)^2/norm(zeros(nodes,1)-x_star)^2;
        err = norm(x0-optx_100_cons)^2/norm(zeros(nodes,1)-optx_100_cons)^2;
        err_value = [err_value;err];
        mu_value = [mu_value,mu0];
        x_value = [x_value,x0];
        fprintf('%d-th iter\n',k);

    end
end
toc
x_value = x_value';
for i = 1:nodes   
    f_value(:,i) = a(i).*x_value(:,i).^2 - b(i).*x_value(:,i);
end
for i = 1:size(f_value,1)
    sum_f(i,:) = sum(f_value(i,:));
end
% filename1 = 'x_value.xlsx';
% filename2 = 'mu_value.xlsx';
% writematrix(x_value,filename1);
% writematrix(mu_value,filename2);