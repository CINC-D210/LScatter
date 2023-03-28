% function bit_seq=tagdecode(p_est,X_change,X,idx)
%#ok<*NASGU> 
%% 测试数据
% 修正后频域栅格上的1200个值  （其余舍弃不用）（是为了接近的目标）
X_change=yRec_tag(:,3);     

% 原本频域栅格上的1200个值  
X=txGrid(:,3);

%% 频域值：从栅格上1200个值 转化为实际2048个IFFT模块 输入的值
N=2048;len=1200;

tmp = complex(zeros(2048,1));
tmp(N/2-len/2+1:N/2, :, :) = X(1:len/2, :, :);  % N=2048  ;len=1200       % 425-1024
tmp(N/2+2:N/2+1+len/2, :, :) = X(len/2+1:len, :, :);    %1025为直流 填0    %  1026-1625        
tmp = [tmp(N/2+1:N, :, :); tmp(1:N/2, :, :)];   %前一半和后一半交换  ？？？
X_in=tmp;  %X_in为原本IFFT前的2048个点（时，频栅格上）是交换了的

tmp = complex(zeros(2048,1));
tmp(N/2-len/2+1:N/2, :, :) = X_change(1:len/2, :, :);  % N=2048  ;len=1200       % 425-1024
tmp(N/2+2:N/2+1+len/2, :, :) = X_change(len/2+1:len, :, :);    %1025为直流 填0    %  1026-1625        
tmp = [tmp(N/2+1:N, :, :); tmp(1:N/2, :, :)];   %前一半和后一半交换  ？？？
X_tag=tmp;  %X_tag为TAg修正后的2048个点（时，频栅格上）
X_valid=X_tag(idx);
% X_tag中的0是人为添加的，实际不是0，但在接收机的dsp中没有收集那些值）

%% 这个OFDM符号中 原本的时域值：
x=ifft(X_in);    %x = step(dsp.IFFT, X_in);   %x为原本IFFT后生成的2048个时域上的值
x_const=[x(1:p_est-1); x(p_est+1200:2048)];
x0=x(p_est:p_est+1200-1);
% % % X_tag-DFT_matrix*x; %验证DFT矩阵的正确性
%% 解码

%得到 Fourier矩阵：
DFT_matrix = dftmtx(2048);  %FFT_matrix=fft(eye(2048));     % DFT矩阵  秩2048
DFT_matrix_idx=DFT_matrix(idx,:);  %有效频率上的那些矩阵（取出了一DFT矩阵的一些行）  （选择行）  秩1200
DTT_matrix_idx_const=DFT_matrix_idx(:,[1:p_est-1,p_est+1200:2048]);  %  有效矩阵中*的是1的那些值   （选择列）秩532
%rank(DTT_matrix_idx_const)=532
%rank(real(DTT_matrix_idx_const))=511
%rank(imag(DTT_matrix_idx_const))=510

W=DTT_matrix_idx_const*x_const;  

Y=X_valid-W;  %Y=F_idx x  的Y（减去常量的1200个频域值）
F_idx=DFT_matrix_idx(:, p_est:p_est+1200-1);   %DFT矩阵中 选出子载波行，去除与自变量无关的列的 子阵  秩739
%F_idx 秩739

% % % %%%验证了 Y=F_idx*x（x为x_real，也就是原本的1200时域值x_weizhi分别乘+1或者-1的新序列）
% % %%%%  基本可以说明和准确的 +-1序列差不多的序列  也会很大程度改变误差  (仅仅改变1位 +-1   0.5749->19.4160 )
% % code=BIT(:,1);
% % % %code(2)=-code(2);
% % x_real=x0.*code;  %BIT 的第一列是第三个OFDM符号上调至的TAG信息
% % error=F_idx*(x_real)-Y;  %  1200个子载波上各自的误差
% % sum(abs(error))  %取误差的模  求和


%% 方程组Ax=b的系数

A=F_idx;  %rank(A)=739
b=Y;
x0; %为原x值
% det(A);  

code=BIT(:,1);   %BIT 的第一列是第三个OFDM符号上调至的TAG信息
x_real=x0.*code;  %x_real 为正确的x值

%%  求解方程组

%1 转化为  By=b  %求+-1 向量 y
B=A*diag(x0);  %rank(B)=737
%验证想法正确
aaa=B*code ;
b;

%2 转化为Cx=d;  %x为0/1的向量
C=-2*B;%复矩阵秩 737     实部秩：1074  虚部秩：1073
d=b-B*ones(1200,1);
% rank(C)=737
% rank(real(C))=1074
% rank(imag(C))=1073

%伪逆
tic
aaa=pinv(real(C))*real(d);
aaa_decode=ones(1200,1);
aaa_decode(find(aaa<0.5))=0; %从虚部解的方程
toc

%验证：
code_01=code;
code_01(find(code==1))=0;
code_01(find(code==-1))=1;
C*code_01;
code_01*code_01';
%% 
clc
tic
bit_real=ones(1200,1);
% [x_real,flag,relres,iter] = lsqr(real(C), real(d),1e-2); %#ok<*ASGLU> 
x_real = lsqminnorm(real(C), real(d),1e-4,'warn');
bit_real(find(x_real<0.5))=0;  %从实部解的方程
toc

x_imag = lsqminnorm(imag(C), imag(d),1e-4,'warn');


bit_imag=ones(1200,1);           
bit_imag(find(x_real<0.5))=0; %从虚部解的方程

bit_errors_real=length(find (bit_real~=code_01))
bit_errors_imag=length(find (bit_imag~=code_01))

err_idx_real=find (bit_real~=code_01);
err_idx_imag=find (bit_imag~=code_01);

err_idx_aaa=find (aaa_decode~=code_01);
%
error_true=norm((C*code_01-d),2);    %真实解的误差
error_solution=norm((C*bit_imag-d),2);  %所求解的误差


%%   

% rank(A)
% rank(real(A))

rank(B)
rank(real(B))

rank(C)
rank(real(C))











%遍历寻找最优解  时间不可受
% tic
% error_min=1e3;
% for cnt=1:2^1200-1
% a=dec2bin(8) ;
% b = double(a)-'0';
% tmp=[zeros(1,1200-length(b)), b];
% tmp(find(tmp==1))=-1;
% tmp(find(tmp==0))=1;
% x_real=x_weizhi.*tmp'; 
% aaaa=Y-F_idx*(x_real);
% error=sum(abs(aaaa));
% if(error<error_min)
%     error_min=error;
%     a_best=a;
% end
% end
% toc
% a

% end