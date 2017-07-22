% Requires Matlab function quadprog

% Description:
% This functions returns the closest matrix U 
% which is reversible to A_input according to 
% the probability density dist and the Frobenius Norm.
% If your not sure what that means, simply read
% page 1 from "Computing the nearest reversible Markov chain"
% from Adam Nielsen and Marcus Weber.

function U=getClosestSparse(A_input, dist, weighted)

% In order to compute the closest sparse matrix, we need to solve
% the convex optimization problem 0.5 x' Q x + f x  und Cx <=0 1*x = 1
% We solve this by using the build in function quadprog from Matlab
% This code only computs the input for that function, i.e Q, f and C.

% Check if input is valid
[a,n] = size(A_input);
if a~=n
    error('Please insert squere matrix.')
end

if n~= size(dist)
   error('Stationary distritubtion has wrong size') 
end

% weighted is optional
if nargin < 3
   weighted = 0
end

% all entries of dist must be non-zero, otherwise it does not induce a
% scalar product!
if all(dist) ~=1 && weighted ==1
	error('Dist has zero entries - not allowed for reweighting scheme')
end

% compute number of basis vectors (see Proposition 2.1)
tempB = sum(dist(:)==0);
m = (n-1)*n/2 +1 + (tempB-1)*tempB/2

% myBasis is an array of m matricies. It contains the basis vectors of the
% supspace U (see Proposition 2.1)
myBasis = cell(1,m);

index=1;
for r=1:(n-1)
    for s=(r+1):n
        if dist(s)==0 && dist(r)==0
            B=speye(n);            
            B(r,r) = 0;
            B(r,s) = 1;
            myBasis{index}=B;
            index=index+1;
            B(r,r) = 1;
            B(r,s) = 0;
            B(s,s) = 0;
            B(s,r) = 1;
            myBasis{index}=B;
        else
            B=speye(n);
            B(r,s) = dist(s);
            B(s,r) = dist(r);
            B(r,r) = 1 - dist(s);
            B(s,s) = 1 - dist(r);
            myBasis{index}=B;
        end        
        index = index + 1;
    end
end
% Last basis vector is the identiy matrix
myBasis{index}=speye(n);

% Compute D and D^(-1) if weighted scheme
if weighted ~= 1
    weighted = 0;
else
   D  = sparse(diag(dist));
   Di = inv(D);
end


f = zeros(m,1);
Q = zeros(m,m);

if weighted == 0
    % Step 1 from 3: Compute f from convex optimization problem
    for i=1:m
        B=myBasis{i};
        f(i) = -2 *trace(B'*A_input);
    end
    
    % Step 2 from 3: Compute Q from convex optimization problem
    for ii=1:m
        B=myBasis{ii};    
        for jj=1:m
            H=myBasis{jj};   
            %Q(i,j)=2*trace(B'*H);
            t = 2*full( B(:).'*H(:) ); % equivalent to trace(B'*H)!
            Q(ii,jj) = t;
            Q(jj,ii) = t;
        end
    end
else
    % Step 1 from 3: Compute f from convex optimization problem
    for i=1:m
        B=myBasis{i};
        f(i) = -2 *trace(D*B*Di*A_input');
    end
    
    % Step 2 from 3: Compute Q from convex optimization problem
    for ii=1:m
        B=myBasis{ii}; 
        Z = D * B *Di;   
        for jj=1:m
            H=myBasis{jj};   
            %Q(i,j)=2*trace(Z'*H);
            t = 2*full( H(:).'*Z(:) ); % equivalent to trace(H'*Z)!
            Q(ii,jj) = t;
            Q(jj,ii) = t;
        end
    end
end



% Step 3 from 3: Compute C from convex optimization problem
C=-eye(m-1+n,m);
C(m,m)=0;

% We only need to compute rows from m to n+m-1 (n rows).
% Each row is equal to -g_i(j) (see page 5 from article).
for i=1:n
        index=1;
        % iterate through basis v_j , j=1,...,m
        % j is stored as alias index.
        for r=1:(n-1)
            for s=(r+1):n
                if dist(s)==0 && dist(r)==0
                    % iterate through 2 basis vectors
                    % delta^[r,s] and delta^[s,r]
                    if r~=i
                        % Case: else
                        C(m-1+i,index)=-1;
                    else
                        % Case: v_j = delta^[i,s]
                        C(m-1+i,index)=0;
                    end
                    index=index+1;
                    if s~=i
                        % Case: else
                        C(m-1+i,index)=-1;
                    else
                        % Case: v_j = delta^[i,s]
                        C(m-1+i,index)=0;
                    end            
                % Case: v_j = A^[r,i]
                elseif s==i
                    C(m-1+i,index)=-1+dist(r);
                % Case: v_j = A^[i,s]
                elseif r==i 
                    C(m-1+i,index)=-1+dist(s);
                % Case: else, this is r~=i && s~=i
                else
                    C(m-1+i,index)=-1
                end
                index=index+1;
            end
        end
        % v_m = Id, thus always g_i(m)=1  
        C(m-1+i,m)=-1;
end

% Start quadprog
zeroVector = zeros(m-1+n,1);
oneVector  = ones(1,m);

options = optimset('Algorithm','interior-point-convex','MaxIter',200);
sol = quadprog(Q,f,C,zeroVector,oneVector,1,[],[],[],options);  

% Compute U
U=zeros(n,n);
for i=1:m
    U =U+ sol(i)*myBasis{i};
end




