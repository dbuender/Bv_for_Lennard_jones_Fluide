function [outputArg] = My(T_red,x)
%MY Summary of this function goes here
%   Detailed explanation goes here

outputArg = exp(-4/T_red*(x^-12-x^-6))-1;
end

