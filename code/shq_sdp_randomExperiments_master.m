close all
clear
clc


dims = 10:10:50;
% dims = 50;
rank = [3,5,7,10];
num_trials = 100;
% num_trials = 1;

%% HPPCA
results_hppca = zeros(length(dims),length(rank),num_trials);

n = [50,200];
v = [1;4];

fprintf('\n HPPCA')
for i=1:length(dims)
    d = dims(i);
    for j=1:length(rank)
        k = rank(j);
        lambda = linspace(1,4,k);
        
        for t = 1:num_trials
            if(t==1 || mod(t,10)==0)
                fprintf('\n dim: %i, rank: %i, trial: %i',d,k,t);
            end
            U = orth(randn(d,k));
            [M,~] = hppca_problem(U,lambda,n,v);
            [proj_err,Xi_err,cvx_optval,Uhat,X,nu,Z,Y] = solve_sdp_CVX(M);
            results_hppca(i,j,t) = Xi_err;
        end
    end
end

save('results_hppca.mat', 'results_hppca');




%% Sym PSD

results_psdM = zeros(length(dims),length(rank),num_trials);
fprintf('\n Symmetric PSD')
for i=1:length(dims)
    d = dims(i);
    for j=1:length(rank)
        k = rank(j);
        for t = 1:num_trials
            if(t==1 || mod(t,10)==0)
                fprintf('\n dim: %i, rank: %i, trial: %i',d,k,t);
            end
            M = generateRand_PSD_M(d,k,k);
            [proj_err,Xi_err,cvx_optval,Uhat,X,nu,Z,Y] = solve_sdp_CVX(M);
            results_psdM(i,j,t) = Xi_err;
        end
    end
end

save('results_psdM.mat', 'results_psdM');


%% Sums of Brocketts
L = [2,5,10];
results_Brockett = zeros(length(dims),length(rank),length(L),num_trials);
fprintf('\n Sums of Brocketts')
for i=1:length(dims)
    d = dims(i);
    for j=1:length(rank)
        k = rank(j);  
        for l=1:length(L)
            for t=1:num_trials
                nGroups = L(l);
                if(t==1 || mod(t,10)==0)
                    fprintf('\n dim: %i, rank: %i, L: %i, trial: %i',d,k,nGroups,t);
                end
                M = generateRandSumsBrocketts(d,k,k,nGroups);
                [proj_err,Xi_err,cvx_optval,Uhat,X,nu,Z,Y] = solve_sdp_CVX(M);
                results_Brockett(i,j,l,t) = Xi_err;
            end
        end
    end
end

save('results_Brockett.mat', 'results_Brockett');


%% Rand M
results_randM = zeros(length(dims),length(rank),num_trials);
fprintf('\n Random Symmetric')
for i=1:length(dims)
    d = dims(i);
    for j=1:length(rank)
        k = rank(j);
        for t=1:num_trials
            if(t==1 || mod(t,10)==0)
                fprintf('\n dim: %i, rank: %i, trial: %i',d,k,t);
            end
            M = generateRand_M(d,k);
            [proj_err,Xi_err,cvx_optval,Uhat,X,nu,Z,Y] = solve_sdp_CVX(M);
            results_randM(i,j,t) = Xi_err;
        end
    end
end

save('results_randM.mat', 'results_randM');

%% Functions to generate random problem instances


function M = generateRand_PSD_M(d,k,rank)
    
    M = cell(k,1);
    for i=1:k
         A = randn(d,rank);
         M{i} = A*A';
    end
    
end

function M = generateRand_M(d,k)
    
    M = cell(k,1);
    for i=1:k
         A = randn(d);
         M{i} = sym(A);
    end
    
end

function M = generateRandSumsBrocketts(d,k,rank,L)
    
    M = cell(k,1);
    W = rand(k,L);

    for i=1:k
        M{i} = zeros(d);
        for l=1:L
            A = randn(d,rank);
            M{i} = M{i} + W(i,l)*(A*A');
        end
        
    end
    
end

function result = sym(A)
result = 0.5*(A + A');
end