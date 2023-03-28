function [out, Eq]  = Equalizer(in, hD, nVar, EqMode)
%#codegen
switch EqMode
    case 1,
       Eq = ( conj(hD))./((conj(hD).*hD));            % Zero forcing
    case 2,
        Eq = ( conj(hD))./((conj(hD).*hD)+nVar);  % MMSE
    otherwise,
        error('Two equalization mode vaible: Zero forcing or MMSE');
end
out=in.*Eq;