function decod_ifft=tagdecode2048_ifft(p_est,X_valid,X)
%#ok<*NASGU> 
%% 测试数据
%TAG调制后，频域栅格上的2048个值  
% X_valid=yRec_tag2048(:,3);     
% 
% % 原本频域栅格上的1200个值  
% X=txGrid(:,3);

%% 频域值：从栅格上1200个值 转化为实际2048个IFFT模块 输入的值
N=2048;len=1200;

tmp = complex(zeros(N,1));
tmp(N/2-len/2+1:N/2, :, :) = X(1:len/2, :, :);  % N=2048  ;len=1200       % 425-1024
tmp(N/2+2:N/2+1+len/2, :, :) = X(len/2+1:len, :, :);    %1025为直流 填0    %  1026-1625        
tmp = [tmp(N/2+1:N, :, :); tmp(1:N/2, :, :)];   %前一半和后一半交换  ？？？
X_in=tmp;  %X_in为原本LTE时频栅格映射后，  IFFT前的2048个点

%% 这个OFDM符号中 原本的时域值：
x=ifft(X_in);    %x = step(dsp.IFFT, X_in);   %x为原本IFFT后生成的2048个时域上的值
x_const=[x(1:p_est-1); ;x(p_est+len:N)];

x0=x(p_est:p_est+len-1);
% % % X_tag-DFT_matrix*x; %验证DFT矩阵的正确性

%% test IFFT

x_change=ifft((X_valid));
x_change1200=x_change(p_est:p_est+len-1);
x1200=x(p_est:p_est+len-1);

code_diff=(x_change1200-x1200)./x_change1200;
%  code_diff=1-(x1200)./x_change1200;
code_diff_abs=abs(code_diff);

decod_ifft=ones(len,1);
decod_ifft(find(code_diff_abs<1))=0;


% code=BIT(:,1);   %BIT 的第一列是第三个OFDM符号上调至的TAG信息
% x_real=x0.*code;  %x_real 为正确的x值
% code_01=code;
% code_01(find(code==1))=0;
% code_01(find(code==-1))=1;
% 
% idx_error=find (decod_ifft~=code_01);
% length(idx_error) %误码率



end