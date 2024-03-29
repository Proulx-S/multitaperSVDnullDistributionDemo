close all
clear all
%% Load example data
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!! Jump to line 70 to visualize precomputed results !!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
load('exampleData.mat')



%% Compute psd
% Included for the sake of completeness
[psd.vec,psd.f] = fastMtPSD(TP.psd.tp,funTs.vec,paramPsd.Fs,[],nBloc,avTapers);
psd.T = size(funTs.vec,1).*tr;


%% Compute multivariate coherence
% Included for the sake of completeness
cohF = [0 inf]; % frequency range over which to compute coherence
nMode = inf;    % number of modes to keep

[svd.u,...      % left-side (voxels) singular vectors
    svd.s,...   % singular spectrum
    svd.v,...   % right-side (tapers) singular vectors
    svd.coh,... % coherence
    svd.f]...   % frequency vector
    = fastKleinMtSVD(TP.svd.tp,funTs.vec,paramSvd.Fs,funTs.t,cohF,nMode);

svd.T = size(funTs.t,1).*tr;
svd.info = 'taper/vox/mode x mode x freq'; % note on matrix dimensions



%% Estimate null distribution of multivariate coherence
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!! This is where the magic happens !!!!
%!!!! -The permutation begins in fastKleinMtSVD.m,line 46.
%!!!! -The matrix that enters the svd is precomputed at line 52.
%!!!! -All possible permutations are precomputed at line 53.
%!!!! -At line 63, a randomly choosen permutation is applied independently
%!!!!  for each voxel, and the shuffled matrix that results is fed to svd,
%!!!!  all in a single line for speed.
nShuf = 2^7;        % number of row-wise random permutation of columns of the voxel x taper matrix that enters the svd
cohFperm = [0 0.35]; % frequency range over which to compute coherence
if nShuf<100; warning('nShuf<100, width of the confidence interval will be underestimated'); end
memFlag = 0; % if 1, stores only the extremes of the null distribution to save on memory, at the expanse of longer (x4) compute time
[svd.uNullMean,...      % Mean of the null distribution of left-side (voxels) singular vectors
    svd.sNullMean,...   % Mean of the null distribution of singular spectrum
    svd.vNullMean,...   % Mean of the null distribution of right-side (tapers) singular vectors
    svd.cohNullMean,... % Mean of the null distribution of coherence
    svd.fNull,...       % frequency vector of the null distribution
    svd.uNull90,...     % 5th and 95th percentile of the null distribution of left-side (voxels) singular vectors
    svd.sNull90,...     % 5th and 95th percentile of the null distribution of singular spectrum
    svd.vNull90,...     % 5th and 95th percentile of the null distribution of right-side (tapers) singular vectors
    svd.cohNull90]...   % 5th and 95th percentile of the null distribution of coherence
    = fastKleinMtSVD(TP.svd.tp,funTs.vec,paramSvd.Fs,funTs.t,cohFperm,nMode,nShuf,[],memFlag);
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



svd.nPerm = nShuf;
svd.infoTaperPerm = 'taper/vox/mode x mode x freq x perm'; % note on matrix dimensions


% %% Save precomputed data
% save(fullfile(wd,'exampleDataPrecomputed.mat'),'svd','-v7.3')
%
% %% Load precomputed data
% load('exampleDataPrecomputed.mat','svd')


%% Plot things
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!! This is where some more magic happens !!!!
%!!!! Figure1, top panel: Coherence spectra, computed using only the first
%!!!! mode in the numerator (1-mode coherence), or with the first n-mode in
%!!!! the numerator (2-mode to 6-mode). Note that each is normalized to a 0
%!!!! to 1 range.
%!!!! Figure1, bottom panel: Rankness of the data, estimated as ("integer")
%!!!! the number of modes before the first one that does not exceed the
%!!!! 95th percentile of the empirically estimate of the null distribution,
%!!!! and as ("crossing point interpolation") the point where the finely
%!!!! interpolated singular spectrum crosses the 95th percentile. For the
%!!!! latter, whenever the first singluar value is already below the 95th
%!!!! percentile, rankness is set to 0.
%!!!! Figure 2 to N-1: Singluar spectra at frequencies indicated by green
%!!!! vertical lines.
%!!!! Figure N, top panel: "1st-mode" coherence as in Figure1, along with
%!!!! ("low-rank") coherence including the first n significant modes in the
%!!!! numerator, n being frequency-specific, and ("low-rank (upsampled)")
%!!!! computed in a similar way but with n estimated on finely interpolated
%!!!! singular spectra. Note that in "low-rank (upsampled)", rank is set to
%!!!! 1 whenever the first mode does not exceed the 95th percentile.
f0List = [0.048 0.271 0.371 0.937 1.118];
plotPerm4(svd,f0List)
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
