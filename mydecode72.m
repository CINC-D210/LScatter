function decod_bit = mydecode72(p_est,yRec_tag, txGrid1,code_len)
N=128;len=72;
repeat=72/code_len;

tmp = complex(zeros(N,1));
tmp(N/2-len/2+1:N/2, :, :) = txGrid1(1:len/2, :, :);  % N=2048  ;len=1200       % 425-1024
tmp(N/2+2:N/2+1+len/2, :, :) = txGrid1(len/2+1:len, :, :);    %1025为直流 填0    %  1026-1625        
tmp = [tmp(N/2+1:N, :, :); tmp(1:N/2, :, :)];   %前一半和后一半交换  ？？？
ofdm_tx_freq=tmp;  %X_in为原本LTE时频栅格映射后，  IFFT前的2048个点
ofdm_tx_time = ifft(ofdm_tx_freq);
tmp = complex(zeros(N,1));
tmp(N/2-len/2+1:N/2, :, :) = yRec_tag(1:len/2, :, :);  % N=2048  ;len=1200       % 425-1024
tmp(N/2+2:N/2+1+len/2, :, :) = yRec_tag(len/2+1:len, :, :);    %1025为直流 填0    %  1026-1625        
tmp = [tmp(N/2+1:N, :, :); tmp(1:N/2, :, :)];   %前一半和后一半交换  ？？？
ofdm_rx_freq=tmp;  %X_in为原本LTE时频栅格映射后，  IFFT前的2048个点
ofdm_rx_time = ifft(ofdm_rx_freq);

decod_bit = ofdm_rx_time./ofdm_tx_time;
decod_bit = decod_bit(p_est:p_est+code_len*repeat-1);
direct = real(mean(decod_bit));
% decod_bit(find(real(decod_bit)>0))=0;
% decod_bit(find(real(decod_bit)<0))=1;
for index=0:code_len-1
    if(sum(real(decod_bit(index*repeat+1:index*repeat+repeat))))>direct
        decod_bit(index*repeat+1:index*repeat+repeat)=0;
    else
        decod_bit(index*repeat+1:index*repeat+repeat)=1;
    end
end
end
    