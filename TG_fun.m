function [fun] = TG_fun(num)
%传统发电机组 cost function
a = [1.42 1.94 4.82]*1e-3;
b = [7.2 7.85 7.97];
c = [510 310 78];

fun = @(x) a(num)*x^2 + b(num)*x + c(num);
end

