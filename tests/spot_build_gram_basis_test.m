function testPass = spot_build_gram_basis_test

    testPass = [];
    
    x = msspoly('x',2);
    u = msspoly('u',1);

    f = u(1)*x(1)^4+x(1)^2*x(2)^2-u(1)*x(2)^2;

    [F,g,pow] = DecompPoly(f,u(1));

    m = spot_build_gram_basis(pow,[],F,g);

    if any(m~=[1,1])
        testPass = [testPass,0];
        spotTestFail('Test failed.')
    else
        testPass = [testPass,1]; 
    end
    
    
    f = u(1)*x(1)^4+x(1)^2-u(1)*x(2)^2;
    [F,g,pow] = DecompPoly(f,u(1));
    m = spot_build_gram_basis(pow,[],F,g);
    
    if any(m~=[1,0])
        testPass = [testPass,0];
        spotTestFail('Test failed.')
    else
        testPass = [testPass,1]; 
    end
    
    
    u = msspoly('u',2);
    f = u(1)*x(1)^4+x(1)^2-u(1)*x(2)^2+u(2)*x(1)^6+u(2)*x(2)^3;
    [F,g,pow] = DecompPoly(f,u);
    m = spot_build_gram_basis(pow,[],F,g);
   
    if any(m~=[1,0])
        testPass = [testPass,0];
        spotTestFail('Test failed.')
    else
        testPass = [testPass,1]; 
    end
    
    
    u = msspoly('u',2);
    f = u(2)*x(2)^3;
    [F,g,pow] = DecompPoly(f,u);
    m = spot_build_gram_basis(pow,[],F,g);
    if ~isempty(m)
        testPass = [testPass,0];
        spotTestFail('Test failed.')
    else
        testPass = [testPass,1]; 
    end 
       
    u = msspoly('u',2);
    x = msspoly('x',3);
    
    f = (u(1)-u(2))*x(1)^4 + u(1)*x(2)^4-u(1)*x(1)^4*x(2)^2+u(2)*x(3)^4+x(3)^2;
    [F,g,pow] = DecompPoly(f,u);
    m = spot_build_gram_basis(pow,[],F,g);
    if any(m~=[0,0,1])
        testPass = [testPass,0];
        spotTestFail('Test failed.')
    else
        testPass = [testPass,1]; 
    end 
 
    
end


function  spotTestFail(string)
    error(string)
end


