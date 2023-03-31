function p2 = prmsDLSCH(cRate,maxIter, fullDecode, p)
p2.trellis = poly2trellis(4, [13 15], 13);
if (cRate >= 1) || (cRate <= 0)
    error('Wrong coding rate');
end
p2.cRate = cRate;
modType = 0.5*p.Qm;
TBLenVec = zeros(1, length(p.numPDSCHBits)); 
C = zeros(1, length(TBLenVec));  Kplus = zeros(1, length(C));
for i = 1:length(TBLenVec)
    TBLenVec(i) = getTBsizeRMC(modType, p2.cRate, p.Nrb, ...
                             p.numLayPerCW, p.numPDSCHBits(i));
    [C(i), ~, Kplus(i)] = lteCblkSegParams(TBLenVec(i));
end
p2.TBLenVec = TBLenVec;
p2.maxTBLen = max(p2.TBLenVec);
p2.maxC = max(C);
p2.minC = min(C);
p2.maxKplus = max(Kplus);
p2.maxIter = maxIter;
p2.fullDecode = fullDecode;