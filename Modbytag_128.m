function [out,offset_int,offset_phi,preamble_p,BIT,BIT01] = Modbytag_128(in)
%输入：时域信号（1ms,14个OFDM符号）(2个符号作为估计信道，p的开销)（估计信道时用全1（））
%输出: 被TAG调制过的时域信号，实际的p,估信道的preamble

% in=txSig;
nFft=128;
nSc=72;
nCp=9;
nCpL=10;

preamble_p=randi([0,1],nSc,1); %估计p的preamble
N_symbol=14*1;                    %本函数中，TAG调制的OFDM符号的个数（在这么多的OFDM符号中认为p不变，信道不变）
N_sample=length(in);         %(2048*15)*1;  %本函数中in的点的个数
% N_sample=2048*15*1;

err_max=nFft-nSc+1-nCp ;  %设计上容忍的最大误差 （防止覆盖后面的CP重复部分，解决多径）
err_min=3;  %硬件所能达到的最小误差   (指的是CP之后的整数p)
offset_int=randi([err_min err_max])  
%延迟的整数倍basic-unit time（CP之后，值得的是有用OFDM符号内的偏移数）
offset_phi=rand(1)*2*pi % 延迟的相位

% offset_int=6;  %测试test
% offset_phi=pi*0.8;  %测试test

tag_bit_01= randi([0,1],nSc,N_symbol);
% tag_bit_01=zeros(nSc,N_symbol); %test
tag_bit_01(:,1)=zeros(nSc,1);  %第一个OFDM符号内的1200比特全部嵌入0  用作估计1200信道
tag_bit_01(:,2)=preamble_p;  %第2个OFDM符号内的1200比特preamble_p作为导频， 估计2048信道及偏移p

tag_ejf=tag_bit_01;  %比特  0 1
tag_ejf(find(tag_ejf==1))=-1;  %将比特1映射成-1（反相位）；
tag_ejf(find(tag_ejf==0))=1;  %将比特0映射成1（同相不变））；
% tag_f=tag_ejf.*exp(j*offset_phi); %加上相偏
tag_f=tag_ejf; 

tag_all=ones(nFft+nCp,N_symbol);

tag_all(nCp+offset_int:nCp+offset_int+nSc-1,:)=tag_f;  %(2048+144)*14  
%先都按照常規CP  但实际上14个符号中的第1，8个符号应该是160长的CP  故后面再加16个0 
%实际硬件中也要考虑这部分延迟
tag_all_serial=tag_all(:);     %转串行  30688
tag_all_serial=[ones(nCpL-nCp,1);tag_all_serial(1:end/2);ones(nCpL-nCp,1);tag_all_serial(end/2+1:end)];


%s-f-tag  TAG调制后信号
out=in.*tag_all_serial.*exp(j*offset_phi);  %TAG反射后的信号  %加上相偏


BIT=tag_ejf(:,3:end);
BIT01=tag_bit_01(:,3:end);
end
