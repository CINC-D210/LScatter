function Mtx=f_repeatMtx(all,nRepeat)

% 返会一个重复矩阵，总共all行，重复nRepeat一变
%总共all行， all/nRepeat列，没列都有nRepeat个1
vec10=zeros(all,1);
vec10(1:nRepeat)=ones(nRepeat,1);
Mtx=zeros(all,all/nRepeat);
for tmp=1:all/nRepeat  
     Mtx(:,tmp)=circshift(vec10,nRepeat*(tmp-1));
end

end