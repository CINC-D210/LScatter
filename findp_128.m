function  [p_est,H128_est]= findp_128(X_change,preamble_p,X,X_tag_rec)
%输入：修正后的频域值 前导码  实际的频域值

% %% 测试数据
% %%修正后恢复的频域栅格上的1200个值  （其余舍弃不用）（是为了接近的目标）
% X_change=yRec_tag(:,2);     

%%原本频域栅格上的1200个值    （填0成2048点频域后，映射成2048时域）
% X=txGrid(:,2);     %这里应该是*了增益的频率值 （Rx时已将Tx的增益抵消了，这里没有增益就是频域上的值）
%过了信道的频域值   1ms内第2个OFDM符号
% X_tag_rec=rxGrid_128(:,2);

%% 从LTE栅格上1200个值 转化为实际2048个IFFT模块 输入的值
N=128;len=72;nCp=9;
tmp = complex(zeros(N,1));
tmp(N/2-len/2+1:N/2, :, :) = X(1:len/2, :, :);  % N=2048  ;len=1200       % 425-1024
tmp(N/2+2:N/2+1+len/2, :, :) = X(len/2+1:len, :, :);    %1025为直流 填0    %  1026-1625        
tmp = [tmp(N/2+1:N, :, :); tmp(1:N/2, :, :)];   %前一半和后一半交换  ？？？
X_in=tmp;  %X_in为LTE原本的  时频栅格映射到IFFT前的128个点
% idx=ExpungeFrom(1:N,find(X_in==0));

tmp = complex(zeros(N,1));
tmp(N/2-len/2+1:N/2, :, :) = X_change(1:len/2, :, :);  % N=2048  ;len=1200       % 425-1024
tmp(N/2+2:N/2+1+len/2, :, :) = X_change(len/2+1:len, :, :);    %1025为直流 填0    %  1026-1625        
tmp = [tmp(N/2+1:N, :, :); tmp(1:N/2, :, :)];   %前一半和后一半交换  ？？？
X_tag=tmp;  %X_tag为TAg初步修正后的128个点（时，频栅格上）
% X_tag中的0是人为添加的，实际不是0，但在接收机的dsp中没有收集那些值）

x=ifft(X_in);    %x = step(dsp.IFFT, X_in);   %x为原本LTE  IFFT后生成的2048个时域上的值

%% 找p
seq=preamble_p;
seq(find(seq==1))=-1;  %将比特1映射成-1（反相位）；
seq(find(seq==0))=1;  %将比特0映射成1（同相不变））；
s_128=ones(N,1);
s_128(1:len)=seq;

S=[];
for tmp=0:N-len+1-nCp  %现在是1-48种可能的偏移量  20230307
     S(:,tmp+1)=circshift(s_128,tmp);
end

nump=size(S,2);
x_list_fft = zeros(N,nump);
for index=1:nump 
    x_list_fft(:,index)=multiAndFilter(x,S(:,index));
end
% x_list=repmat(x,1,nump).*S;
% x_list_fft=fft(x_list);

X_tag_rp=repmat(X_tag,1,nump);
compare_index = [2:37 93:128];
error_vec=abs(x_list_fft(compare_index,:)-X_tag_rp(compare_index,:));  %这里误差=|实际值-真实值|  

error_sum=sum(error_vec,1);%将误差累加
[~,p_est]=min(error_sum); 

%% 求各个子载波上信道值
X_tag_theory = multiAndFilter(x,S(:,p_est));
%过了信道的频域值  X_tag_rec是 1ms内第2个OFDM符号
H128_est=X_tag_rec./X_tag_theory; %得到了128个子载波上的信道值
% figure
% plot(abs(H128_est))

% %% 验证时域符号生成的正确性   说明了时域上发送的数据txSig，比IFFT产生的时域数据大2048/sqrt(1200)倍
% test=x.*(2048/sqrt(1200));  %时域上实际发送的 功率增大的信号
% yuan=txSig(2048+160+144+1:2048+160+144+1+2047);  %验证时域符号生成的正确性 
% mean(yuan-test)


end

