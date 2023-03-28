clc;clear;close

N=128
x=randn([N,1])+1j*randn([N,1]);
fai=pi/4;
x_phi=x.*exp(1j*fai);%时域上每个点都*e_jfai

fft1=fft(x);
fft2=fft(x_phi);

t=fft2./fft1;
mean(t)
exp(1j*fai)

%% 说明时域上