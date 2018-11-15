function vOut = trimData(vFormants,vBW,nLowerBound,nBW)

nn = 1;
vOut = [];

for kk = 1:length(vFormants)
   if ((vFormants(kk) > nLowerBound) && (vBW(kk) < nBW))
       vOut(nn) = vFormants(kk);
       nn = nn+1;
   end
end

end