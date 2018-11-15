function vOut = guaranteeLength(vIn,nLenDesired)

% Guarantee that input vector vIn has length nDesired

vIn = vIn(:);

nLen = numel(vIn);
if nLen > nLenDesired
    vOut = vIn(1:nLenDesired);
elseif (nLen == nLenDesired)
    vOut = vIn;
else
    vOut = [vIn;zeros(nLenDesired-nLen,1)];
end

end