function bit_seq=tagdecode128(p_est,X_valid,X)
%#ok<*NASGU> 
%% 测试数据
% %TAG调制后，频域栅格上的2048个值  
% X_valid=yRec_tag2048(:,3);     
% 
% % 原本频域栅格上的1200个值  
% X=txGrid(:,3);

%% 频域值：从栅格上72个值 转化为实际128个IFFT模块 输入的值
N=128;len=72;

tmp = complex(zeros(N,1));
tmp(N/2-len/2+1:N/2, :, :) = X(1:len/2, :, :);  % N=2048  ;len=1200       % 425-1024
tmp(N/2+2:N/2+1+len/2, :, :) = X(len/2+1:len, :, :);    %1025为直流 填0    %  1026-1625        
tmp = [tmp(N/2+1:N, :, :); tmp(1:N/2, :, :)];   %前一半和后一半交换  ？？？
X_in=tmp;  %X_in为原本LTE时频栅格映射后，  IFFT前的2048个点

%% 这个OFDM符号中 原本的时域值：
x=ifft(X_in);    %x = step(dsp.IFFT, X_in);   %x为原本IFFT后生成的2048个时域上的值
x_const=[x(1:p_est-1); x(p_est+len:N)];

x0=x(p_est:p_est+len-1);
% % % X_tag-DFT_matrix*x; %验证DFT矩阵的正确性


%% 解码

%得到 Fourier矩阵：
DFT_matrix = dftmtx(N);  %FFT_matrix=fft(eye(2048));     % DFT矩阵  秩2048
DTT_matrix_const=DFT_matrix(:,[1:p_est-1,p_est+len:N]);  %  有效矩阵中*的是1的那些值   （选择列）秩532
%rank(DTT_matrix_idx_const)=532
%rank(real(DTT_matrix_idx_const))=511
%rank(imag(DTT_matrix_idx_const))=510

W=DTT_matrix_const*x_const;  
Y=X_valid-W;  %Y=F_idx x  的Y（减去常量的1200个频域值）
F_idx=DFT_matrix(:, p_est:p_est+len-1);   %DFT矩阵中 选出子载波行，去除与自变量无关的列的 子阵  秩739

%% 方程组Ax=b的系数

A=F_idx;  
b=Y;
% x0; %为原x值
% 
% code=BIT(:,1);   %BIT 的第一列是第三个OFDM符号上调至的TAG信息
% x_real=x0.*code;  %x_real 为正确的x值

%%  求解方程组

%1 转化为  By=b  %求+-1 向量 y
B=A*diag(x0);  %rank(B)=737
% % % %验证想法正确
% % % aaa=B*code ;
% % % b;

%2 转化为Cx=d;  %x为0/1的向量
C=-2*B;%复矩阵秩 737     实部秩：1074  虚部秩：1073
d=b-B*ones(len,1);
% C = C([2:37,93:128],:);
% d = d([2:37,93:128],:);

%伪逆求

% tic%%用实部
% decode=inv(real(C)'*real(C))*real(C)'*(real(d));
% toc

   %复数矩阵计算
   warning('off')
decode=inv(C'*C)*C'*(d);
% decode_01=ones(len,1);
% decode_01(find(abs(decode)<0.5))=0; %从解的方程


bit_seq=ones(len,1);
bit_seq(find(abs(decode)<0.5))=0; %从解的方程


% % %验证：
% % code_01=code;
% % code_01(find(code==1))=0;
% % code_01(find(code==-1))=1;
% % err_idx_real=find (aaa_decode~=code_01);
% % 
% % bit_imag=code;
% % bit_imag(find(code==1))=0;
% % bit_imag(find(code==-1))=1;

% % error_true=norm((C*code_01-d),2);    %真实解的误差
% % error_solution=norm((C*bit_imag-d),2);  %所求解的误差

end