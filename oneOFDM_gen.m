function  x= oneOFDM_gen(in)

%输入为频域1200个值 生成其时域值（现在未加CP） 2022.12.1
% in=txGrid(:,2);
persistent hIFFT;
if isempty(hIFFT)
    hIFFT = dsp.IFFT;
end
[len, numSymb, numLayers] = size(in);
% N assumes 15KHz subcarrier spacing
N = 2048;
cpLen0 = 144;

tmp = complex(zeros(N, numSymb, numLayers));
% Pack data, add DC, and reorder
tmp(N/2-len/2+1:N/2, :, :) = in(1:len/2, :, :);  % N=2048  ;len=1024        % 513-1024
tmp(N/2+2:N/2+1+len/2, :, :) = in(len/2+1:len, :, :);    %1025为直流 填0    %  1026-1537        
tmp = [tmp(N/2+1:N, :, :); tmp(1:N/2, :, :)];   %前一半和后一半交换  ？？？
% IFFT processing
x = step(hIFFT, tmp);
x = x.*(N/sqrt(len));   %功率归一化

end