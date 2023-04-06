function decode_01=decode72(p_est,X_change,X,idx,nRepeat)
%只利用72子载波的数据，解调72/nRepeat个TAG比特

%% 测试数据
% 修正后频域栅格上的1200个值  （其余舍弃不用）（是为了接近的目标）
% X_change=yRec_tag(:,3);     
% 
% % 原本频域栅格上的1200个值  
% X=txGrid(:,3);

%% 频域值：从栅格上1200个值 转化为实际2048个IFFT模块 输入的值
N=128;len=72;

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
x_const=[x(1:p_est-1); x(p_est+len:N)];
x0=x(p_est:p_est+len-1);

%% 解码

%得到 Fourier矩阵：
DFT_matrix = dftmtx(N);  %FFT_matrix=fft(eye(2048));     % DFT矩阵  秩2048
DFT_matrix_idx=DFT_matrix(idx,:);  %有效频率上的那些矩阵（取出了一DFT矩阵的一些行）  （选择行）  秩1200
DTT_matrix_idx_const=DFT_matrix_idx(:,[1:p_est-1,p_est+len:N]);  %  有效矩阵中*的是1的那些值   （选择列）秩532


W=DTT_matrix_idx_const*x_const;  
Y=X_valid-W;  %Y=F_idx x  的Y（减去常量的1200个频域值）

F_idx=DFT_matrix_idx(:, p_est:p_est+len-1);   %DFT矩阵中 选出子载波行，去除与自变量无关的列的 子阵  秩739

%% 方程组Ax=b的系数
A=F_idx;  %rank(A)=739
b=Y;
% x0; %为原x值
% rank(A)

% code=BIT(:,1);   %BIT 的第一列是第三个OFDM符号上调至的TAG信息
% x_real=x0.*code;  %x_real 为正确的x值
%%  求解方程组
%1 转化为  By=b  %求+-1 向量 y
B=A*diag(x0);  %rank(B)=737
%验证想法正确

% testB=B*code-b;

%2 转化为Cx=d;  %x为0/1的向量
C=-2*B;%复矩阵秩 737     实部秩：1074  虚部秩：1073
d=b-B*ones(len,1);

% testC=C*code_01-d;


%% 考虑重复的情况：
%3 转化为Dx=d 
%D: 72*(72/nRepeat)  x:(72/nRepeat)*1    e:72*1

vec10=zeros(72,1);
vec10(1:nRepeat)=ones(nRepeat,1);
repeatMtx=zeros(72,72/nRepeat);
for tmp=1:72/nRepeat  
     repeatMtx(:,tmp)=circshift(vec10,nRepeat*(tmp-1));
end

D=C*repeatMtx;
d;

% tic
decode=inv(D'*D)*D'*(d);
decode_01=ones(72/nRepeat,1);
decode_01(find(abs(decode)<0.5))=0; %从解的方程
% toc

% err_idx=find (decode_01~=bit01(:,1));
% bit_error=length(err_idx);

% rank(C)
% rank(real(C))


end