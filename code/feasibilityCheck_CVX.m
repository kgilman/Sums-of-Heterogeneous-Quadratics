function [result,nu] = feasibilityCheck_CVX(U,M)

%%% Inputs:
%  U: d x k orthonormal matrix that is a critical point on the Stiefel
%  manifold
%  M: k-length cell of dxd symmetric PSD matrices
%%% Outputs:
%  returns 0 if U is certified to be globally optimal; returns 1 otherwise
%  returns feasible nu for dual certificate; otherwise returns array of
%  NaN's.


%%% Check that U is a critical point on the manifold
[d,k] = size(U);
tol = 1e-7;
AU = zeros(d,k);
for i=1:k
    AU(:,i) = M{i}*U(:,i);
end
if(norm(AU,2) > tol)
    result = 1;
    nu = NaN*ones(k,1);
end

Ik = eye(k);
Lam = zeros(k);
for i=1:k
    ei = Ik(:,i);
    Lam = Lam + (U'*M{i}*U*ei)*ei';
end

Lambda = sym(Lam);

if(norm(Lambda - Lam,'fro') / norm(Lambda,'fro') > tol || ~ispsd(Lambda,tol))
    result = 1;
    nu = NaN*ones(k,1);
end


%%% run the feasibility problem in CVX
    cvx_begin sdp quiet
        cvx_precision best
%         cvx_precision default
    %     cvx_solver mosek
        variable nu(k,1)

    %     maximize(sum(nu));

        subject to
        for i=1:k
           sym(U*(Lambda - diag(nu))*U') + nu(i)*eye(d) - sym(M{i}) >= 0;
        end

        Lambda - diag(nu) >= 0;

    cvx_end

    if(strcmp(cvx_status,'Solved') || strcmp(cvx_status, 'Inaccurate/Solved'))
        result = 0;
    elseif(strcmp(cvx_status, 'Inaccurate/Solved'))
        result = checkDualFeasibleSolns(U,M,nu);
    else
        result = 1;
    end
end



function result = sym(A)
    result = 0.5*(A + A');
end
