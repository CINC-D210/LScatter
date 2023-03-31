function y = OFDMTx(in, prmLTE)  
%生成1子帧（1ms）的时域信号：把1200个子载波，14个OFDM符号内的频域值 映射成时域信号
% in为资源网格：
% 第一维度：子载波1200 维度的复数（列）
% 第二维度：OFDM符号数 
% 第三维度：层数
% prmLTE.N;  %IFFT点数
% prmLTE.cpLen0;    %第一个CP长度 160 
% prmLTE.cpLenR;    %余下6个CP长度 144

% in=txGrid;
% prmLTE=prmLTEPDSCH;

%#codegen
persistent hIFFT;
if isempty(hIFFT)
    hIFFT = dsp.IFFT;
end
[len, numSymb, numLayers] = size(in);
% N assumes 15KHz subcarrier spacing
N = prmLTE.N;
cpLen0 = prmLTE.cpLen0;
cpLenR = prmLTE.cpLenR;
slotLen = (N*7 + cpLen0 + cpLenR*6);
subframeLen = slotLen*2;  
tmp = complex(zeros(N, numSymb, numLayers));
% Pack data, add DC, and reorder
tmp(N/2-len/2+1:N/2, :, :) = in(1:len/2, :, :);  % N=2048  ;len=1200       % 425-1024
tmp(N/2+2:N/2+1+len/2, :, :) = in(len/2+1:len, :, :);    %1025为直流 填0    %  1026-1625        
tmp = [tmp(N/2+1:N, :, :); tmp(1:N/2, :, :)];   %前一半和后一半交换  ？？？
% IFFT processing
x = step(hIFFT, tmp); %每                                                                                                                                                                                                                                                                                                                一列的两边为0，中间是频域（交换后的频域值）
x = x.*(N/sqrt(len));   %**********功率归一化*******
% Add cyclic prefix per OFDM symbol per antenna port 
% and serialize over the subframe (equal to 2 slots)
% For a subframe of data
y = complex(zeros(subframeLen, numLayers));
for j = 1:2 % Over the two slots
    % First OFDM symbol
    y((j-1)*slotLen+(1:cpLen0), :) = x((N-cpLen0+1):N, (j-1)*7+1, :);
    y((j-1)*slotLen+cpLen0+(1:N), :) = x(1:N, (j-1)*7+1, :);

    % Next 6 OFDM symbols
    for k = 1:6
        y((j-1)*slotLen+cpLen0+k*N+(k-1)*cpLenR+(1:cpLenR), :) = x(N-cpLenR+1:N, (j-1)*7+k+1, :);
        y((j-1)*slotLen+cpLen0+k*N+k*cpLenR+(1:N), :) = x(1:N, (j-1)*7+k+1, :);
    end
end