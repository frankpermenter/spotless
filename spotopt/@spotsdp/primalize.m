function [prout,G,h] = primalize(prin)
%
%  [prout,G,h] = primalize(prin)
%
%  Converts a program into the standard primal with free
%  variables form via the introduction of slack variables.
%
%  The matrices G,h are constructed so that
%
%  prg.variables = G*spPrg.variables + h.
%
    prout = prin;
    prout.posCnst = [];
    prout.lorCstr = [];
    prout.lorCstrDim = [];
    prout.rlorCstr = [];
    prout.rlorCstrDim = [];
    prout.psdCstr = [];
    prout.psdCstrDim = [];
    
    [pos,dim] = prin.posConstraints();
    if ~isempty(dim) && dim ~= 0
        [prout,slack] = prout.newPos(dim);
        prout = prout.withEqs(pos-slack);
    end
    
    [lor,dim] = prin.lorConstraints();
    if ~isempty(dim)
        [prout,slack] = prout.newLors(dim);
        prout = prout.withEqs(lor-slack);
    end
    
    [rlor,dim] = prin.rlorConstraints();
    if ~isempty(dim)
        [prout,slack] = prout.newRLors(dim);
        prout = prout.withEqs(rlor-slack);
    end
    
    [psd,dim] = prin.psdConstraints();
    if ~isempty(dim)
        [prout,slack] = prout.newPSDs(dim);
        prout = prout.withEqs(psd-slack);
    end
    
    h = zeros(size(prin.variables));
    
    [var,pow,Coeff] = decomp(prout.variables);
    mtch = match(var,prin.variables);
    G = Coeff(:,mtch)';
end