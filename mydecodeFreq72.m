function decod_bit=mydecodeFreq72(p_est,X_recovered,X_tx,code_len)
%% 频域值：从栅格上72个值 转化为实际128个IFFT模块 输入的值
N=128;len=72;
repeat = len/code_len;
tmp = complex(zeros(N,1));
tmp(N/2-len/2+1:N/2, :, :) = X_tx(1:len/2, :, :);  % N=2048  ;len=1200       % 425-1024
tmp(N/2+2:N/2+1+len/2, :, :) = X_tx(len/2+1:len, :, :);    %1025为直流 填0    %  1026-1625        
tmp = [tmp(N/2+1:N, :, :); tmp(1:N/2, :, :)];   %前一半和后一半交换  ？？？
X_in=tmp;  %X_in为原本LTE时频栅格映射后，  IFFT前的2048个点

X_recovered = [X_recovered(len/2+1:len);X_recovered(1:len/2)];

%% 这个OFDM符号中 原本的时域值：
x=ifft(X_in);    %x = step(dsp.IFFT, X_in);   %x为原本IFFT后生成的2048个时域上的值
x_const=[x(1:p_est-1); x(p_est+len:N)];

x0=x(p_est:p_est+len-1);
% % % X_tag-DFT_matrix*x; %验证DFT矩阵的正确性

%% 解码

%得到 Fourier矩阵：
DFT_matrix = dftmtx(N);  %FFT_matrix=fft(eye(2048));     % DFT矩阵  秩2048
DTT_matrix_const=DFT_matrix(:,[1:p_est-1,p_est+len:N]);  %  有效矩阵中*的是1的那些值   （选择列）秩532

W=DTT_matrix_const*x_const;  
Y=X_recovered-W([2:len/2+1,end-len/2+1:end]);  %Y=F_idx x  的Y（减去常量的72个频域值）
F_idx=DFT_matrix([2:len/2+1,N-len/2+1:N], p_est:p_est+len-1);   %DFT矩阵中 选出子载波行，去除与自变量无关的列的 子阵  秩739

%% 方程组Ax=b
A=F_idx*diag(x0);  
b=Y;
%% C*x=d x取0/1
% C=-2*A;
% d=b-A*ones(len,1);
% decode=inv(C'*C)*C'*(d);
% decod_bit=ones(len,1);
% decod_bit(find(abs(decode)<0.5))=0; %从解的方程
decod_bit=lsqminnorm(A,b);
% decod_bit=inv(A)*b;
for index=0:code_len-1
    if(sum(real(decod_bit(index*repeat+1:index*repeat+repeat))))>0
        decod_bit(index*repeat+1:index*repeat+repeat)=0;
    else
        decod_bit(index*repeat+1:index*repeat+repeat)=1;
    end
end

end