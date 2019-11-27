
function [mykspace_zf, NCol, NLin] = zerofillkspace(mykspace,NCol0,NLin0, k_0fill)
    
if (k_0fill > 0)
     
    m     = nextpow2(NLin0) + k_0fill;
    n     = nextpow2(NLin0) + k_0fill;
    padx  = ((2^m) - NLin0)/2;
    pady  = ((2^n) - NLin0)/2;

    mykspace_zf = padarray(mykspace,[padx,pady],0); 
    NCol = size(mykspace_zf,1);  
    NLin = size(mykspace_zf,2);
    w = window2(NCol,NLin,@hamming);
    
    mykspace_zf = mykspace_zf.* w;
    
else 
    NCol = NCol0; 
    NLin = NLin0;
    mykspace_zf = mykspace;
    
end

end