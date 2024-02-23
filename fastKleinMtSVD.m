function [u,s,v,c,f,u5to95,s5to95,v5to95,c5to95] = fastKleinMtSVD(tapers,data,Fs,tvec,cohFpass,nMode,nPerm,verbose,memFlag)
% tapers must be [tapers x 1 x time]
if ~exist('verbose','var') || isempty(verbose); verbose = 1; end
if ~exist('nPerm','var'); nPerm = 0; end
if ~exist('memFlag','var'); memFlag = false; end
[NC,C]=size(data); % size of data
[NK,K]=size(tapers,[3 1]); % size of tapers
if NK~=NC; error('length of tapers is incompatible with length of data'); end;
NT = NC; clear NC NK
if ~exist('nMode','var') || isempty(nMode)
    nMode = 1;
elseif nMode == inf
    nMode = min([C K]);
end
if ~exist('cohFpass','var') || isempty(cohFpass)
    cohFpass = [0 0.2];
end

%% Time vector
if ~exist('tvec','var') || isempty(tvec)
    tvec=1/Fs *(0:NT-1)';
elseif isduration(tvec)
    tvec = seconds(tvec);
end
tvec=squeeze(tvec*2*pi*1i);

%% Frequencies
nfft=max(2^(nextpow2(NT)+0),NT);
freqs=getfgrid(Fs,nfft,[0 Fs/2]);
NF=length(freqs);
freqs = permute(freqs,[1 3 2 4]);
if cohFpass(2)==inf; cohFpass(2) = freqs(end); end
fInd = freqs>=cohFpass(1) & freqs<=cohFpass(2);
f = freqs(fInd);
nf = nnz(fInd);

%% Perform SVD (vox X taper) at each frequencies
if ~nPerm
    %%%%%%%%%%%%%%%%%
    %% Compute SVD %%
    %%%%%%%%%%%%%%%%%
    [u,s,v] = pagesvd(reshape(permute(data,[2 1])*reshape(permute(tapers,[3 1 2]).*exp(-freqs(fInd).*tvec),[NT K*nf]),[C K nf]),"econ","vector");
    c = s.^2./sum(s.^2,1); % coherence
    u5to95 = []; s5to95 = []; v5to95 = []; c5to95 = [];
else 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Estimate null distribution of SVD results %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Prepare data matrix (mat) for fast permutations
    % mat = reshape(permute(data,[2 1])*reshape(permute(tapers,[3 1 2]).*exp(-freqs(fInd).*tvec),[NT K*nf]),[C K nf]);
    mat = permute(reshape(permute(data,[2 1])*reshape(permute(tapers,[3 1 2]).*exp(-freqs(fInd).*tvec),[NT K*nf]),[C K nf]),[3 2 1]);
    allPerm = perms(1:K); allPermN = size(allPerm,1);
    
    if ~memFlag
        %% Store the full distribution and summarize after (faster but more memory)
        tic
        %%% Generate full distribution
        sTmp =         zeros(min([C K]),1,nf,nPerm,'single') ;
        uTmp = complex(zeros(C,min([C K]),nf,nPerm,'single'));
        vTmp = complex(zeros(K,min([C K]),nf,nPerm,'single'));
        for permInd = 1:nPerm
            [uTmp(:,:,:,permInd),sTmp(:,:,:,permInd),vTmp(:,:,:,permInd)] = pagesvd(permute( reshape( mat(:,allPerm(randi(allPermN,C,1),:)' + (0:K:C*K-1)) , [nf K C]) , [3 2 1] ),"econ","vector");
        end
        %%% Compute first-mode coherence
        cTmp = sTmp.^2./sum(sTmp.^2,1); % [mode X 1 X freq]

        %%% Summarize ditribution
        s = mean(sTmp,4); u = mean(uTmp,4); v = mean(vTmp,4); c = mean(cTmp,4);
        s5to95 = cat(4,prctile(abs(sTmp),5,4),prctile(abs(sTmp),95,4)); clear sTmp
        u5to95 = cat(4,prctile(abs(uTmp),5,4),prctile(abs(uTmp),95,4)); clear uTmp
        v5to95 = cat(4,prctile(abs(vTmp),5,4),prctile(abs(vTmp),95,4)); clear vTmp
        c5to95 = cat(4,prctile(abs(cTmp),5,4),prctile(abs(cTmp),95,4)); clear cTmp
        toc
    else
        %% Store partial distribution (less memory but longer)
        tic
        u = zeros(C,         K,nf,'single');
        s = zeros(K,         1,nf,'single');
        v = zeros(min([C K]),K,nf,'single');
        c = zeros(K,         1,nf,'single');
        nPermExtr = ceil(nPerm*0.05); if nPermExtr<2; nPermExtr = 2; end
        u95 = zeros(C,         K,nf,nPermExtr,'single');
        s95 = zeros(K,         1,nf,nPermExtr,'single');
        v95 = zeros(min([C K]),K,nf,nPermExtr,'single');
        c95 = zeros(K,         1,nf,nPermExtr,'single');
        u5 = repmat(realmax('single'),[C,         K,nf,nPermExtr]);
        s5 = repmat(realmax('single'),[K,         1,nf,nPermExtr]);
        v5 = repmat(realmax('single'),[min([C K]),K,nf,nPermExtr]);
        c5 = repmat(realmax('single'),[K,         1,nf,nPermExtr]);
        % figure('WindowStyle','docked');
        for permInd = 1:nPerm
            [uTmp,sTmp,vTmp] = pagesvd(permute( reshape( mat(:,allPerm(randi(allPermN,C,1),:)' + (0:K:C*K-1)) , [nf K C]) , [3 2 1] ),"econ","vector");
            cTmp = sTmp.^2./sum(sTmp.^2,1);
            u = u+uTmp; s = s+sTmp; v = v+vTmp; c = c+cTmp;
            u95 = replacePerm(u95,uTmp,'top'); s95 = replacePerm(s95,sTmp,'top'); v95 = replacePerm(v95,vTmp,'top'); c95 = replacePerm(c95,cTmp,'top');
            u5 = replacePerm(u5,uTmp,'bot'); s5 = replacePerm(s5,sTmp,'bot'); v5 = replacePerm(v5,vTmp,'bot'); c5 = replacePerm(c5,cTmp,'bot');
            % histogram([squeeze(s5(1,1,round(end/2),:)); squeeze(s95(1,1,round(end/2),:))]); drawnow
        end
        u = u./nPerm; s = s./nPerm; v = v./nPerm; c = c./nPerm;
        %%% Summarize ditribution
        u5to95 = cat(4,max(u5,[],4),min(u95,[],4)); clear u5 u95
        s5to95 = cat(4,max(s5,[],4),min(s95,[],4)); clear s5 s95
        v5to95 = cat(4,max(v5,[],4),min(v95,[],4)); clear v5 v95
        c5to95 = cat(4,max(c5,[],4),min(c95,[],4)); clear c5 c95
        toc
    end
end
% u: spatial singular vector [vox X mode X frequency]
% s: singular values [mode X 1 X frequency]
% v: taper singular vector [taper X mode X frequency]
% c: coherence [mode X 1 X frequency]

%% Delete higher modes
if nMode<min([C K])
    u(:,nMode+1:end,:,:) = [];
    v(:,nMode+1:end,:,:) = [];
    if exist('u5to95','var'); u5to95(:,nMode+1:end,:,:) = []; end
    if exist('v5to95','var'); v5to95(:,nMode+1:end,:,:) = []; end
end

