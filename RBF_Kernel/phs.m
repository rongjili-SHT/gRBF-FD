classdef phs
    
    methods(Static)
        function v = phi(r,kappa)
            v = r.^(2*kappa+1);
        end
        
        %%% Beltrami 2nd order derivate
        function dd = beltrami(r,kappa,dim)
            dd = ((2*kappa+1)*(2*kappa - 1) + dim*(2*kappa+1)) .* r.^(2*kappa-1) ;
        end
    end
end

