close all
clear all

%% Load example data
wd = pwd;
if exist(fullfile(wd,'exampleData.mat'),'file')
    % From this folder
    load(fullfile(wd,'exampleData.mat'))
elseif exist(fullfile([wd '_data'],'exampleData.mat'),'file')
    % From another "_data" this folder
    wd = [wd '_data'];
    load(fullfile(wd,'exampleData.mat'))
end


%% Compute psd
[psd.vec,psd.f] = fastMtPSD(TP.psd.tp,funTs.vec,paramPsd.Fs,[],nBloc,avTapers);
psd.T = size(funTs.vec,1).*tr;


%% Compute multivariate coherence
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

svd.nPerm = nShuf;
svd.infoTaperPerm = 'taper/vox/mode x mode x freq x perm'; % note on matrix dimensions

save exampleDataPrecomputed svd -v7.3

f0List = [0.048 0.271 0.371 0.937 1.118];
plotPerm4(svd,f0List)
