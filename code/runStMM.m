function [Uhat,fxn_log] = runStMM(M,U0,params,fxn)

H.U = U0;

[d,k] = size(U0);

% F0 = F(H.U,M);
fxn_log = zeros(params.niters + 1,1);

fxn_log(1) = fxn(H.U);
for iter=1:params.niters
    dFdU = gradF(H.U,M);

    fo_err = norm((dFdU - H.U*(dFdU'*H.U)),'fro');
    if(fo_err <= params.break_tol)
        fxn_log = fxn_log(1:iter);
        break
    end
    
    [U,~,V] = svd(dFdU,0);
    H.U = U*V';
    if(norm(H.U' * H.U - eye(k)) > 1e-6)
        error("Loss of orthonormality in StGA")
    end

    fxn_log(iter + 1) = fxn(H.U);

    
end

Uhat = H.U;

end


function grad = gradF(U,M)
    grad = zeros(size(U));
    for k=1:size(U,2)
        grad(:,k) = M{k}*U(:,k);
    end
end

