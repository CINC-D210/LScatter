function hD=gridResponse_interpolate(hp, Nrb, Nrb_sc, Ndl_symb)
% Average over the two same Freq subcarriers, and then interpolate between
% them - get all estimates and then repeat over all columns (symbols).
% The interpolation assmues NCellID = 0.
% Time average two pilots over the slots, then interpolate (F)
% between the 4 averaged values, repeat for all symbols in sframe
hD = complex(zeros(Nrb*Nrb_sc, Ndl_symb*2));
N=size(hp,2);
Separation=6;
Edges=[0,5;3,2;0,5;3,2];
Symbol=[1,5,8,12];
% First: Compute channel response over all resource elements of OFDM symbols 0,4,7,11
for n=1:N
    Edge=Edges(n,:);
    y = InterpolateCsr(hp(:,n),  Separation, Edge);
    hD(:,Symbol(n))=y;
end
% Second: Interpolate between OFDM symbols {0,4} {4,7}, {7, 11}, {11, 13}
for m=[2, 3, 4, 6, 7]
    alpha=0.25*(m-1);
    beta=1-alpha;
    hD(:,m)    = beta*hD(:,1) + alpha*hD(:,  5);
    hD(:,m+7) =beta*hD(:,8) + alpha*hD(:,12);
end
