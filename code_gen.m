clc; clear
code_len = 24; %%只调4*18个bit
repeat=3;
offset = (72-code_len*repeat)/2;

tag_preamble_ori = repmat([1 0]',code_len/2,1);%%0是0，1是pi
tag_preamble=repmat(tag_preamble_ori,1,repeat).';%方波调制信号
tag_preamble=tag_preamble(:);

tag_data_ori = randi([0,1],code_len,10);
tag_data = zeros(code_len*repeat,10);
for index=1:10
    tag_data_repeat_temp=repmat(tag_data_ori(:,index),1,repeat).';%方波调制信号
    tag_data(:,index)=tag_data_repeat_temp(:);
end

cp = 9+1;
all_data = zeros(1920,1);
for index = 4:14 %12 unused, 34 for phi, p&channel
    if index==8
        cp=cp+1;
    end
    if index == 4
        all_data((index-1)*137+cp+1+offset:(index-1)*137+cp+code_len*repeat+offset) = tag_preamble;
    else
        all_data((index-1)*137+cp+1+offset:(index-1)*137+cp+code_len*repeat+offset) = tag_data(:,index-4);
    end
end
all_data_char = char(all_data'+'0')
save('tag_code.mat',"tag_preamble","tag_data","all_data_char");

% tag_preamble=preamble_p;
% tag_data = Bit01;
% save('tag_code.mat',"tag_preamble","tag_data");

% W = zeros(128*2,1);
% w = vertcat(zeros((128-code_len)/2,1),tag_preamble,zeros((128-code_len)/2,1));
% for i = 0:127
%     if w(i+1) == 0
%         phase = 1;
%     else 
%         phase = -1;
%     end
%     W(i*2+1) = w(i+1);
%     W(i*2+2) = -w(i+1);
% end
% plot(abs(fft(W)))

% tag_preamble = load("tag_code.mat").tag_preamble;
% plot(abs(fft(randi([-1,1],2048))));
% tag_data = load("tag_code.mat").tag_data;
% 
% all_data = zeros(1920,1);
% all_data(137+cp+1:137+cp+code_len) = tag_preamble;
% for index = 3:14
%     if index>=8
%         cp=cp+1;
%     end
%     all_data((index-1)*137+cp+1:(index-1)*137+cp+code_len)=tag_data(:,index-2);
% end
% all_data_char = char(all_data'+'0')
% save('tag_code.mat',"tag_preamble","tag_data","all_data_char");
% % var(tag_data)
% % sum(tag_data,1)
