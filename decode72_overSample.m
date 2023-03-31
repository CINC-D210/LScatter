function decode_01=decode72_overSample(p_est,X_change,X,idx,nRepeat,nOvers)
%只利用72子载波的数据，解调72/nRepeat个TAG比特

% %% 测试数据
% %TAG调制后频域栅格上的1200个值  （其余舍弃不用）（是为了接近的目标）
% X_change=yRec_tag(:,14);     
% 
% % 原本频域栅格上的1200个值  
% X=txGrid(:,14);
% nRepeat=code_len;
% nOvers=10;

%% 
N=128;len=72;
tmp = complex(zeros(N,1));
tmp(N/2-len/2+1:N/2, :, :) = X(1:len/2, :, :);  % N=2048  ;len=1200       % 425-1024
tmp(N/2+2:N/2+1+len/2, :, :) = X(len/2+1:len, :, :);    %1025为直流 填0    %  1026-1625        
tmp = [tmp(N/2+1:N, :, :); tmp(1:N/2, :, :)];   %前一半和后一半交换  ？？？
X_in=tmp;  %X_in为原本IFFT前的72个点（时，频栅格上）是交换了的

%OverSample后的频域，中间是0
lteFre_O = vertcat(X_in(1:64),zeros((nOvers-1)*128,1),X_in(65:128));

%% 这个OFDM符号中 空中信号的时域值：
x=ifft(lteFre_O)*1;     %x为LTE的时域信号，IFFT后到时域上的值（过采样的） ***系数**
%% 频域值：从栅格上72个值 转化为实际128个IFFT模块 输入的值


%USRP采集到的频域值
tmp = complex(zeros(N,1));
tmp(N/2-len/2+1:N/2, :, :) = X_change(1:len/2, :, :);  % N=2048  ;len=1200       % 425-1024
tmp(N/2+2:N/2+1+len/2, :, :) = X_change(len/2+1:len, :, :);    %1025为直流 填0    %  1026-1625        
tmp = [tmp(N/2+1:N, :, :); tmp(1:N/2, :, :)];   %前一半和后一半交换  ？？？
X_tag=tmp;  %X_tag为TAg修正后的2048个点（时，频栅格上）
X_valid=X_tag(idx);  %%%TAG调制后，USRP收到的频域72值（为了逼近的方程右边的向量（目标））
% X_tag中的0是人为添加的，实际不是0，但在接收机的dsp收集的那些值不准确？？）

%% 解码

%得到 Fourier矩阵：

DFT_matrix = dftmtx(N*nOvers);     % DFT矩阵 
DFT_matrix2=DFT_matrix*diag(x)*f_repeatMtx(N*nOvers,nOvers);  %*diag(x)*重复矩阵

% %%测试
% code=Bit01(:,cnt-4);  %测试
% code_zf1=ones(72,1);
% code_zf1(find(code==1))=-1;  %1 变-1     0（其余）均是1
% code_128=ones(128,1);
% code_128(p_est:p_est+72-1)=code_zf1;
% 
% f_1280=DFT_matrix2*code_128;
% f_128=dftmtx(N)*diag(ifft(X_in))*code_128;     
% 
% X_valid; %TAG调制后的72


DFT_matrix_idx=DFT_matrix2(idx,:);  %有效频率上的那些矩阵（取出了一DFT矩阵的一些行）  （选择行）  
DTT_matrix_idx_const=DFT_matrix_idx(:,[1:p_est-1,p_est+len:N]);  %  有效矩阵中*的是1的那些值   （选择列）

x_const=ones(128-72,1);
W=DTT_matrix_idx_const*x_const;  
Y=X_valid-W;  %Y=F_idx x  的Y（减去常量的1200个频域值）

F_idx=DFT_matrix_idx(:, p_est:p_est+len-1);   %DFT矩阵中 选出子载波行，去除与自变量无关的列的 子阵  秩739

%% 求解方程组  Bx=b
%1 转化为  By=b  %求+-1 向量 y
B=F_idx; 
b=Y;
%验证想法正确

% %测试
% code=Bit01(:,cnt-4);
% code_zf1=ones(72,1);
% code_zf1(find(code==1))=-1;  %1 变-1     0（其余）均是1
% testB=B*code_zf1-b;



%2 转化为Cx=d;  %x为0/1的向量
C=-2*B;%复矩阵秩 737     实部秩：1074  虚部秩：1073
d=b-B*ones(len,1);

% testC=C*code-d;


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

%误码率测试
% Bit_down=downsample(Bit01(:,cnt-4),nRepeat);
% err_idx=find (decode_01~=Bit_down);
% bit_error=length(err_idx)

end
% end