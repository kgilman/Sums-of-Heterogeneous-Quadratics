close all
clear
clc

thresh = 1e-5;


%% HPPCA
load results_hppca.mat results_hppca

% test = max(results_hppca,[],3)

array = 1 - sum(results_hppca > thresh,3) ./ 100;

dims = {'d=10';'d=20';'d=30';'d=40';'d=50'};
r3 = array(:,1);
r5 = array(:,2);
r7 = array(:,3);
r10 = array(:,4);


T = table(r3,r5,r7,r10,'RowNames',dims)
writetable(T,'results_hppca.dat','WriteRowNames',true) 

%% Sym PSD


load results_psdM.mat results_psdM

% test = max(results_hppca,[],3)

array = 1 - sum(results_psdM > thresh,3) ./ 100;

dims = {'d=10';'d=20';'d=30';'d=40';'d=50'};
r3 = array(:,1);
r5 = array(:,2);
r7 = array(:,3);
r10 = array(:,4);


T = table(r3,r5,r7,r10,'RowNames',dims)
writetable(T,'results_psdM.dat','WriteRowNames',true) 

%% randM
load results_randM.mat results_randM

% test = max(results_hppca,[],3)

array = 1 - sum(results_randM > thresh,3) ./ 100;

dims = {'d=10';'d=20';'d=30';'d=40';'d=50'};
r3 = array(:,1);
r5 = array(:,2);
r7 = array(:,3);
r10 = array(:,4);


T = table(r3,r5,r7,r10,'RowNames',dims)
writetable(T,'results_randM.dat','WriteRowNames',true)

%% Sums of Brocketts
load results_Brockett.mat results_Brockett

% test = max(results_hppca,[],3)

array = 1 - sum(results_Brockett > thresh,4) ./ 100;

dims = {'d=10, L=2'; 'd=10, L=5'; 'd=10, L=10';...
        'd=20, L=2'; 'd=20, L=5'; 'd=20, L=10';...
        'd=30, L=2'; 'd=30, L=5'; 'd=30, L=10';...
        'd=40, L=2'; 'd=40, L=5'; 'd=40, L=10';...
        'd=50, L=2'; 'd=50, L=5'; 'd=50, L=10'};
r3 = vec(array(:,1,:));
r5 = vec(array(:,2,:));
r7 = vec(array(:,3,:));
r10 = vec(array(:,4,:));


T = table(r3,r5,r7,r10,'RowNames',dims)
writetable(T,'results_Brockett.dat','WriteRowNames',true)

%% HPPCA Extended

load results_hppca_extended.mat results_hppca

% test = max(results_hppca,[],3)

ns = [5, 10, 20, 50, 100, 200];
% n = [50,200];
vs = [1;2;3;4];
thresh = 1e-5;

p = 4;
for m=1:length(ns)

    array = 1 - sum(squeeze(results_hppca(:,:,m,p,:)) > thresh ,3) ./ 100;

    dims = {'d=10';'d=20';'d=30';'d=40';'d=50'};
    r3 = array(:,1);
    r5 = array(:,2);
    r7 = array(:,3);
    r10 = array(:,4);


    T = table(r3,r5,r7,r10,'RowNames',dims)
    writetable(T,'results_hppca.dat','WriteRowNames',true) 

end

disp("Across v")
m = 2;
for p=1:length(vs)

    array = 1 - sum(squeeze(results_hppca(:,:,m,p,:)) > thresh ,3) ./ 100;

    dims = {'d=10';'d=20';'d=30';'d=40';'d=50'};
    r3 = array(:,1);
    r5 = array(:,2);
    r7 = array(:,3);
    r10 = array(:,4);


    T = table(r3,r5,r7,r10,'RowNames',dims)
    writetable(T,'results_hppca.dat','WriteRowNames',true) 

end
