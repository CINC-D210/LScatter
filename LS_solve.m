%% SDR例程

n=1200;
H=real(B);
y=real(b);
s0=ones(n,1);
L=[H', -H*y; -y'*H, y'*y];
% W=[s0*s0',s0;s0',1];

 cvx_solver Gurobi
%% 
cvx_begin sdp quiet;
variable W(n+1,n+1) hermitian semidefinite;
minimize(real(trace(L*W)));
subject to
diag(W)==1;
W>=0;
cvx_end

[ SOLVER, SOLVER_LIST ] = cvx_solver 