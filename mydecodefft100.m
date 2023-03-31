function decod_bit = mydecodefft100(p_est,yRec_tag, txGrid1,code_len)
N=128;len=72;
data_num=code_len;
repeat=72/code_len;

tmp = complex(zeros(N,1));
tmp(N/2-len/2+1:N/2, :, :) = txGrid1(1:len/2, :, :);  % N=2048  ;len=1200       % 425-1024
tmp(N/2+2:N/2+1+len/2, :, :) = txGrid1(len/2+1:len, :, :);    %1025为直流 填0    %  1026-1625        
tmp = [tmp(N/2+1:N, :, :); tmp(1:N/2, :, :)];   %前一半和后一半交换  ？？？
ofdm_tx_freq=tmp;  %X_in为原本LTE时频栅格映射后，  IFFT前的2048个点
ofdm_tx_time = ifft(ofdm_tx_freq);
% 
% tmp = complex(zeros(N,1));
% tmp(N/2-len/2+1:N/2, :, :) = yRec_tag(1:len/2, :, :);  % N=2048  ;len=1200       % 425-1024
% tmp(N/2+2:N/2+1+len/2, :, :) = yRec_tag(len/2+1:len, :, :);    %1025为直流 填0    %  1026-1625        
% tmp = [tmp(N/2+1:N, :, :); tmp(1:N/2, :, :)];   %前一半和后一半交换  ？？？
% ofdm_rx_freq=tmp;  %X_in为原本LTE时频栅格映射后，  IFFT前的2048个点
% ofdm_rx_time = ifft(ofdm_rx_freq);

freq_72s = zeros(72,2^data_num);
ZEROS1 = zeros(p_est-1,1);
ZEROS2 = zeros(57-p_est,1);
for code_index=0:2^data_num-1
    code_bin = dec2bin(code_index);
    code_24 = ones(data_num,1)*-1;
    code_24(end-length(code_bin)+1:end) = 2*(code_bin-'0')-1;
    code_temp=repmat(code_24,1,repeat).';
    code_72=code_temp(:);
    [~,freq_72s(:,code_index+1)]=multiAndFilter(ofdm_tx_time,[ZEROS1; code_72; ZEROS2]);%构造频域调制信号
end
freq72_rp=repmat(yRec_tag,1,2^data_num);%复制得到接收到的频域调制信号
error_vec=abs(freq_72s-freq72_rp);  %这里误差=|实际值-真实值|  
error_sum=sum(error_vec,1);%将误差累加
[~,code_est]=min(error_sum); 
code_bin = dec2bin(code_est-1);
code_24 = ones(data_num,1)*-1;
code_24(end-length(code_bin)+1:end) = 2*(code_bin-'0')-1;
code_temp=repmat(code_24,1,repeat).';%方波调制信号
decod_bit=code_temp(:);
decod_bit(find(decod_bit==1))=0;
decod_bit(find(decod_bit==-1))=1;
end
    