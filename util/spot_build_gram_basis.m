function exponent_m = spot_build_gram_basis(pow)
    
    
    exponent_p_monoms = pow;
    csclasses={1:size(pow,2)};
    exponent_m = monomialgeneration(exponent_p_monoms,csclasses);
    
    options = sdpsettings;
    options.verbose = 0;
    temp=sdpvar(1,1);
    tempops = options;
    tempops.solver = 'cdd,glpk,*';  % CDD is generally robust on these problems
    tempops.verbose = 0;
    tempops.saveduals = 0;
    [aux1,aux2,aux3,LPmodel] = export(set(temp>0),temp,tempops);  
    %disp('Reducing Monomials.');
    exponent_m = monomialreduction(exponent_m,exponent_p_monoms,options,csclasses,LPmodel);
    exponent_m = exponent_m{1};
end