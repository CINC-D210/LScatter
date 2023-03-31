function v = mydet(A)

[m,n] = size(A);
if m~=n
    error('This is not a square matrix!')          % 如果A不是方阵，报错
else
    if n>1
        v = 0;
        for j=1:n
            M1j = A;
            M1j(1,:) = [];
            M1j(:,j) = [];                         % M1j是划掉A的第1行和第j列剩下的n-1阶方阵
            v = v + (-1)^(1+j)*A(1,j)*mydet(M1j);  % 利用拉普拉斯展开，递归调用mydet
        end
    else
        v = A(1,1);                                % 最终A降为1阶，行列式即A(1,1)元素
    end
end

% % 作者：Boo Radley
% % 链接：https://www.zhihu.com/question/505046568/answer/2269765537
% % 来源：知乎
% % 著作权归作者所有。商业转载请联系作者获得授权，非商业转载请注明出处。