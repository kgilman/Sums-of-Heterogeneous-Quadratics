function [Uhat,fxn_log] = runStGA(M,U0,params,fxn)

H.U = U0;

[d,k] = size(U0);

% F0 = F(H.U,M);
fxn_log = zeros(params.niters + 1,1);
% fxn_log(1) = F(H.U,M);
fxn_log(1) = fxn(H.U);
for iter=1:params.niters
    [H,terminate] = updateU(H,M,params);
    if(norm(H.U' * H.U - eye(k)) > 1e-6)
        error("Loss of orthonormality in StGA")
    end

%     fxn_log(iter+1) = F(H.U,M);
    fxn_log(iter + 1) = fxn(H.U);
%     if(abs(Ft - F0) <= params.tol)
% %         break
%     else
%         F0 = Ft;
%     end

    if(terminate)
        fxn_log = fxn_log(1:iter+1);
        break
    end
end

Uhat = H.U;

end

