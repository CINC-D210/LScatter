function [freq128, freq72] = multiAndFilter(ofdm_time,tag_time)
%multiAndFilter 对两段长128的，USRP和tag时域信号相乘并在频域理想滤波

nRepeat = 100; 

% % %%%test1
% ofdm_time = eNodeBOutput(11:10+128);
% data_num = 10;
% tag_time=randi([0 1],data_num,1);
% tag_time(find(tag_time==0)) = -1;
% tag_time = vertcat(ones((128-data_num)/2,1),tag_time,ones((128-data_num)/2,1));
% % %%%test
% data_num = 20;
% tag_time=randi([0 1],data_num,1);
% tag_time(find(tag_time==0)) = -1;
% tag_time=repmat(tag_time,1,3).';%方波调制信号
% tag_time=tag_time(:);
% tag_time = vertcat(ones((128-data_num*3)/2,1),tag_time,ones((128-data_num*3)/2,1));
% tag_time=ones(128,1);
% tag_time(64)=1;
% %%%test

ofdm_freq = fft(ofdm_time);%%已知时域lte信号，变道freq dom
%理想低通滤波器
ofdm_freq_n_freq = vertcat(ofdm_freq(1:64),zeros((nRepeat-1)*128,1),ofdm_freq(65:128));
ofdm_freq_n_time = ifft(ofdm_freq_n_freq);%空中的usrp时域信号（没过信道）

tag_repeat=repmat(tag_time,1,nRepeat).';%方波调制信号
tag_repeat2=tag_repeat(:);

modu_n_time = ofdm_freq_n_time.*tag_repeat2;
modu_n_freq = fft(modu_n_time);

freq128 = [modu_n_freq(1:64);modu_n_freq(end-63:end)];%理论上TAG修正后频域的值（没过信道）
freq72 = [modu_n_freq(end-35:end); modu_n_freq(2:37)];%理论上TAG修正后频域的值（没过信道）
% figure
% plot(abs(freq2_5))

% power_per = sum(abs(freq2_5).*abs(freq2_5))/sum(abs(modu_n_freq).*abs(modu_n_freq))

end
%  plot(abs(fft(ofdm_time)))
%  plot(abs(freq2_5))