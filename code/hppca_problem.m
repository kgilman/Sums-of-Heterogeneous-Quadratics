function [M,Y] = hppca_problem(U,lambda,n,v)
    [d,k] = size(U);
    L = length(v);
    k = length(lambda);
    U = U(:,1:k);
    LF = U*diag(sqrt(lambda));

    Y = cell(L,1);
    for l=1:L
        Y{l} = LF*randn(k,n(l)) + sqrt(v(l))*randn(d,n(l));
%         Y{l} = LF*randn(k,n(l));
    end

    w = @(l,i) lambda(i) / v(l) / (lambda(i) +v(l));
%     w = @(l,i) 1;

    %%%%%%%%%%%%% HPPCA-generated M_i %%%%%%%%%%%%%%%%%%%%
    M = cell(k,1);
    for i=1:k
        M{i} = zeros(d);
        for l=1:L
            M{i} = M{i} + w(l,i)*Y{l}*Y{l}';         %%hppca
        end
    end
end