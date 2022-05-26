function ans=ispsd(A,tol)
    d = eig(A);
%     eps = 1e-12;
%     tol = length(d)*eps(max(d));
%     tol = 1e-7;
    if(all(d >= -tol))
        ans= 1;
    else
        ans= 0;
    end
end