function [mpow] = spot_build_gram_basis(pow,mpow,F,g)
  
    pow = full(pow);
    
    if nargin < 2 || isempty(mpow)
       mpow = spot_exponent_bound_polytope(pow);  
    else
       mpow = full(mpow); 
    end
       
    numHyperPlanes = 2000;
    [mpow] = RandomPrune(pow,mpow,numHyperPlanes);
    [mpow] = DiagConsistent(pow,mpow);

    if nargin > 2
        mpow = FacialReduction(pow,mpow,F,g);
    end

end
    
     
%Remove monomials m if 2*m is outside convhull(pow) using
%2*numHyperPlanes random seperating hyperplanes
function mpow = RandomPrune(pow,mpow,numHyperPlanes)

    M = size(pow,2);
    
    for i=1:numHyperPlanes
        
        w = randn(M,1);
        thres1 = max(pow*w);
        thres2 = min(pow*w);
 
        y = 2*mpow*w;
        
        %keep mpow if 2*mpow is in convhull(pow)
        mpow = mpow(y <= thres1 & y >= thres2,:);
               
        %Quit if we've pruned mpow enough.
        if (size(mpow,1) < 100)
            break;
        end
        
    end
        
end

%Get the "diagonally" consistent monomials. 2*mpow must be in pow
%or 2*mpow must be cancelled by a cross term
function [mpow] = DiagConsistent(pow,mpow)

    if (isempty(mpow))
       return 
    end

    %sorting and casting seems to speed up intersect()
    dataType = 'int16';
    pow = sortrows(cast(pow,dataType));
    mpow = sortrows(cast(mpow,dataType));
   
    while (1)

        [mext,mextIndx] = ComputeExtremeMonomials(mpow);
        [~,mextDeleteIndx] = setdiff(2*mext,pow,'rows');
        
        if (~isempty(mextDeleteIndx))
            mpow(mextIndx(mextDeleteIndx),:) = [];
        else
            break;
        end

    end

    mpow = double(mpow);
      
end

%mext: monomials that aren't the midpoints of other monomials
%generalized notion of "extreme point"
function [mext,mextIndx,crossTerms] = ComputeExtremeMonomials(mpow)

    if (size(mpow,1) == 1)
       mext = mpow; 
       mextIndx = 1;
       crossTerms = [];
       return; 
    end
    
    crossTerms = CrossTerms(mpow);
    [~,mextIndx] = setdiff(mpow*2,crossTerms,'rows');
    mext = mpow(mextIndx,:);
    
end

%Compute A+A, excluding term A_i + A_j if i = j
function [S] = CrossTerms(A)

    numTerm = size(A,1);
    numVar = size(A,2);
    
    %Get index position of upper triangular part, minus diagonal
    posKeep = find(triu(ones(numTerm,numTerm,'uint8'),1));

    if ~isempty(posKeep)
        %Compute A(i,:)+A(j,:) for all i \ne j
        S=zeros(length(posKeep),numVar,class(A));

        for k=1:numVar
           temp=bsxfun(@plus,A(:,k),A(:,k)');
           S(:,k) = temp(posKeep);
        end
        
    else
        S = [];
    end

end


%Apply algorithm from Permenter, Parrilo CDC 2014.
function [mpow] = FacialReduction(pow,mpow,F,g)

    if ~exist('linprog','file')
        return 
    end

    if (isempty(mpow))
       return 
    end
    
    powinput = pow;
    Finput = F;
    ginput = g;
    
    while(1)
        
        pow = powinput;
        [mext,mextIndx,minkSum] = ComputeExtremeMonomials(mpow);
        
        %The only relevant monomials are given pow and 2*mext.
        %Add 2*mext to pow if necessary and use updated pow to index
        %all computations.
        powadd = setdiff(2*mext,pow,'rows');
        pow = [powinput;powadd];
        F = [Finput;sparse(size(powadd,1),size(Finput,2))];
        g = [ginput;sparse(size(powadd,1),1)];
        

        [~,facePerpIndx] = setdiff(pow,[mext*2;minkSum],'rows');
        [~,faceDualIndx,faceDualMonomIndx] = intersect(pow,mext*2,'rows');

        Ffree =  F(facePerpIndx,:);
        gfree =  g(facePerpIndx,:);
        Fnneg =  F(faceDualIndx,:);
        gnneg =  g(faceDualIndx,:);

        %Build and solve lp
        l = [-Inf*ones(length(facePerpIndx),1);zeros(length(faceDualIndx),1)];
        u = Inf*ones(length(facePerpIndx)+length(faceDualIndx),1);

        sumEq = [zeros(1,length(facePerpIndx)),ones(1,length(faceDualIndx))];

        A = [[Ffree;Fnneg]';[gfree;gnneg]';sumEq];
        b = [zeros(size(A,1)-1,1);1];

        [x,~,flag] = linprog(zeros(length(u),1),[],[],A,b,l,u);

        FLAG_SOLUTION_FOUND = 1;
        eps = 10^-5;
        
        %update or quit
        if (flag == FLAG_SOLUTION_FOUND)
            xdual = x(length(facePerpIndx)+1:end);
            deleteIndx = faceDualMonomIndx(xdual > eps);
            mpow(mextIndx(deleteIndx),:) = [];
            
            if (isempty(mpow))
               break
            end
            
        else
            break;
        end
        
    end
    
    
end









