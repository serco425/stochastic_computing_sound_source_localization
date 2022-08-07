function [OrApproxResult] = MultiOrApproxBinary(vec,LSB)
% Does the Add approximation of OR gate for Unary Bitstreams using Binary
% calculation

OrApproxResult = 0;

if exist('LSB','var') %

    vec = myquant(vec,LSB);
    
   for i = 1:length(vec)
        OrApproxResult = myquant(OrApproxResult + vec(i) - myquant(OrApproxResult*vec(i),LSB),LSB);
   end
else
    
    for i = 1:length(vec)
        OrApproxResult = OrApproxResult + vec(i) - OrApproxResult*vec(i);
    end   
    
end

end