%% DEMO for Interpolation on the manifold of k-GMMs
% Hyunwoo J. Kim, Nagesh Adluru, Monami Banerjee, Baba C. Vemuri, Vikas Singh, 
% Interpolation on the manifold of k component Gaussian Mixture Models (GMMs),
% In International Conference on Computer Vision (ICCV), December 2015. 
%
% Project page.
% http://pages.cs.wisc.edu/~hwkim/projects/k-gmm/
%
% Github repository
%
%
% Github page
%
%

%% Generate random gmm distributions.

clear
rng default
addpath(genpath('.'))
kcomp = 3;
d = 2;
gmm1 = randgmm(d, kcomp);
gmm2 = randgmm(d, kcomp);

%% GMM parameterizations

% L2 normalization can be obtained by 
ngmm1 = l2normalizeGMM(gmm1);
ngmm2 = l2normalizeGMM(gmm2);
fprintf('L2 norm of the normalized GMMs |ngmm1| = %f, |ngmm2| = %f \n', l2normGMM(ngmm1), l2normGMM(ngmm2));

% Square-root reparameteraization can be obtained by taking the square root
% of the PDF values.

%% Distance between GMMs

% L2 distance between GMMs
d_l2 = l2distGMM(gmm1, gmm2);

% L2 distance between normalized GMMs
% d_nl2 = l2distGMM(ngmm1, ngmm2);
d_nl2 = l2ndistGMM(gmm1,gmm2);

% Geodesic distance on the unit Hilbert sphere with the l2-normalized GMMs.
d_ngeo = acos(innerprodGMM(ngmm1, ngmm2));

fprintf('L2 norm of the normalized GMMs d_l2(gmm1,gmm2) = %f, d_nl2(gmm1, gmm2) = %f, d_ngeo(gmm1,gmm2) = %f \n', d_l2, d_nl2, d_ngeo);

%% Gradient descent scheme for $l_2$-distance

clear
close all
gmm1 = obj2structGMM(gmdistribution([-5; 2], ones(1,1,2),[0.7 0.3]));
gmm2 = obj2structGMM(gmdistribution([-4.9 ;2.1;3.1; 4;5], ones(1,1,5),[0.2 0.4 0.1 0.25 0.05]));
gmm3 = obj2structGMM(gmdistribution([-5.1 ;2.5;4;], ones(1,1,3),[0.2 0.4 0.4]));

gmms = {gmm1, gmm2, gmm3};
gmml2mean = l2meanGMMs(gmms);

figure
myplotgmm1D(gmml2mean,[-10,10],[1,0,1],0.01);
title(sprintf('L2 mean of GMMs with k= %d', length(gmml2mean.PComponents)));
maxiter = 100;
k = 3;
tic
[kgmm, fval, gnorm, fval_history, status, gmm0] = gd_gmms_closest_fast(gmml2mean, k, maxiter);
toc
figure
myplotgmm1D(gmm0,[-10,10],[1,0,1],0.01);
title('Random Initialization kGMM');

figure
myplotgmm1D(kgmm,[-10,10],[1,0,1],0.01);
title('Learned GMM');

%% EM algorithm for KL-divergence or cross entory.
% Note that the proposed EM algorithms do not optimize L2 distance.
% But it returns a good solution w.r.t L2 distance too.

option.maxiter = 100;
[gmmbar_em, stats] = em_gmm_closest_full(gmml2mean, k, option)

figure;
hold on
myplotgmm1D(gmmbar_em,[-15 15],  [1,0,1], 0.01);
title('EM GMM fitting.')
hold off

%% Modified EM (our restricted GPMM)
% Note that the proposed EM algorithms do not optimize L2 distance.
% But it returns a good solution w.r.t L2 distance too.

option.getGamma = 'getGamma_l2dist';
option.debug = true;
gmmbar_mem = em_rgpmm(gmml2mean, k, option); 

figure;
hold on
myplotgmm1D(gmmbar_mem,[-15 15],  [1,0,1], 0.01);
title('Modified EM GMM fitting.')
hold off