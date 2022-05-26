function fval = F(U,M)
    fval = 0;
    for k=1:size(U,2)
        fval = fval + U(:,k)'*M{k}*U(:,k);
    end
end