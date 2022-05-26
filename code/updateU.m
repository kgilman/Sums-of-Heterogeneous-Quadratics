
function [H,terminate] = updateU(H,M,params)
%     params = sga.stepsize;
    
%     dFdU = gradF(M.U,M.lam,M.v,Y);
    terminate = 0;
    dFdU = gradF(H.U,M);
    nablaF = dFdU - H.U*(dFdU'*H.U);
    if(norm(nablaF) < params.break_tol)
        terminate = 1;
        return
    end

%     F0 = F(M.U,M.lam,M.v,Y);
    F0 = F(H.U,M);
    FDelta = params.tol * norm(nablaF)^2;
    for m=0:params.maxsearches - 1
%     m = findfirst(IdentityUnitRange(0:params.maxsearches-1)) do m
        Delta = params.contraction^m * params.stepsize;
        if(F(geodesic(H.U,nablaF,Delta),M) >= F0 + Delta * FDelta)
            break;
        end
        if(m == params.maxsearches - 1)
            sprintf("Exceeded maximum line search iterations. Accuracy not guaranteed.")
        end
    end

    Delta = params.contraction^m * params.stepsize;
    H.U = geodesic(H.U,nablaF,Delta);
%     return M;
end

function Unext = geodesic(U,X,t)
    k = size(U,2);

    A = skew(U'*X);
    [Q,R] = qr(X - U*(U'*X));
    Q = Q(:,1:k);
    R = R(1:k,1:k);

    MN = expm(t*[A -R'; R zeros(k,k)]);
    MN = MN(:,1:k);
    M = MN(1:k,:);
    N = MN(k+1:end,:);

    Unext = U*M + Q*N;
end

function result = skew(A)
    result = (A - A')/2;
end

% function grad = gradF(U,lam,v,Y)
%     grad = zeros(size(U));
%     for l=1:length(Y)
%         grad = grad + Y(l) * (Y(l)' * (U * diag(lam./v(l)./(lam + v(l)))));
%     end
% 
% end

function grad = gradF(U,M)
    grad = zeros(size(U));
    for k=1:size(U,2)
        grad(:,k) = M{k}*U(:,k);
    end
end

% function fval = F(U,lam,v,Y)
%     fval = 0;
%     for l=1:length(Y)
%         fval = fval + norm(sqrt(diag(lam./v(l)./(lam + v(l))))*U'*Y(l))^2;
%     end
%     fval = 0.5*fval;
% end

function fval = F(U,M)
    fval = 0;
    for k=1:size(U,2)
        fval = fval + 0.5*U(:,k)'*M{k}*U(:,k);
    end
end

% gradF = @(U,lam,v,Y) sum(Yl * Yl' * U * diag(lam/vl/(lam +vl)) for (Yl,vl) in zip(Y,v))
% F = @(U,lam,v,Y) 1/2*sum(norm(sqrt(diag(lam/vl/(lam + vl)))*U'*Yl)^2 for (Yl,vl) in zip(Y,v))