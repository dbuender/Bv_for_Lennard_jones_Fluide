clc, clear all, close all
N  = 300;
RN = rand(1,N);
RN_n = randn(1,N);
subplot(1,2,1),histogram(RN)
xlabel('\chi'), ylabel('Frecuency')
subplot(1,2,2),histogram(RN_n)
xlabel('\chi'), ylabel('Frecuency')
