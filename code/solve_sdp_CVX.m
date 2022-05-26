function [proj_err,Xi_err,cvx_optval,Uhat,X,nu,Z,Y] = solve_sdp_CVX(M)

k = length(M);
[d,~] = size(M{1});
%%%%%%%%%%%% CVX to solve the SDP %%%%%%%%%%%%%%
cvx_begin sdp quiet
%     cvx_solver mosek
    cvx_precision high
 
    variable X(d,d,k) symmetric
    dual variables nu{k}
    dual variables Z{k}
    dual variable Y
    
    obj = 0;
    for i=1:k
        obj = obj + trace(M{i}*X(:,:,i));
    end
    maximize(obj);
    
    XX = zeros(d,d);
    for i=1:k
        X(:,:,i)  >= 0: Z{i};
        trace(X(:,:,i)) == 1: nu{i};
        XX = XX + X(:,:,i);
    end
    XX <= eye(d): Y;

cvx_end

proj_err = norm(XX - XX*XX,'fro');
Xi_err = 0;
Id = eye(d);
for i=1:k
    Xi_err = Xi_err + 1/k*norm(sort(eig(X(:,:,i)),'descend') - Id(:,1));
end

Uhat = zeros(d,k);
for i=1:k
    [U2,D2] = eig(X(:,:,i));
    [~,idx]=sort(diag(D2),'descend');
    Uhat(:,i) = U2(:,idx(1))*D2(idx(1),idx(1));
end

[Uhat,~] = qr(Uhat,0);

% [Q,~,V] = svd(Uhat,0);

% Uhat = Q*V';


end