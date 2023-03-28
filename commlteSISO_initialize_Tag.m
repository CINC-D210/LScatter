function [prmMdl_Tag] = commlteSISO_initialize_Tag()
%建立结构体变量： PDSCH资源块生成相关  DLSCH解码相关   Mdl信道相关

% Channel parameters
prmMdl_Tag.chanMdl = chanMdl;
prmMdl_Tag.corrLevel = 'frequency-selective-low-mobility';
prmMdl_Tag.chEstOn = 1;  %是否估计信道


end