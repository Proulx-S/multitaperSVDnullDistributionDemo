function [J,f] = fastMtPSD(tapers,data,Fs,Fpass,nBloc,avTapers,verbose)
% tapers must be [1 x 1 x 1 x time x tapers]
sz = size(data,1:4);
NR = sz(4);
NC = sz(1);
C = prod(sz([2 4]));
% [NC,C]=size(data); % size of data
[NK,K]=size(tapers,[4 5]); % size of tapers
if NK~=NC; error('length of tapers is incompatible with length of data'); end;
NT = NC; clear NC NK
if ~exist('Fpass','var') || isempty(Fpass)
    Fpass = [0 inf];
end
if ~exist('nBloc','var') || isempty(nBloc)
    nBloc = 1;
end
if ~exist('avTapers','var') || isempty(avTapers)
    avTapers = false;
    if nBloc>1
        avTapers = true;
    end
end
if ~exist('verbose','var') || isempty(verbose)
    verbose = 2;
end

%% FFT
pad = 0;
NFFT=max(2^(nextpow2(NT)+pad),NT);
if nBloc==1
    % Compute everything at once
    if avTapers
        if verbose>1; disp('averaging across tapers--phase information will be lost'); end
        % Project to taper space -> compute fft -> convert to power -> average
        % across tapers -> sqrt for consistency
        J=sqrt(mean(abs(reshape(fft(reshape(data(:,:).*permute(tapers,[4 1 5 2 3]),[NT C*K]),NFFT)/Fs,[NFFT C K])).^2,3));
    else
        % Project to taper space -> compute fft -> keep data from all tapers
        J=reshape(fft(reshape(data(:,:).*permute(tapers,[4 1 5 2 3]),[NT C*K]),NFFT)/Fs,[NFFT C K]);
    end
else
    % Split computations across nBloc blocs of roughly equal number of
    % voxels, averaging across tapers and deleting data between each bloc
    % to save on memory
    if verbose>1; disp('averaging across tapers--phase information will be lost'); end
    Cdist = repmat(round(C/nBloc),[1 nBloc]); Cdist(end) = Cdist(end) + (C - sum(Cdist));
    data = mat2cell(data(:,:),NT,Cdist);
    J = cell(size(data));
    for bloc = 1:nBloc
        tic
        if verbose; disp(['bloc ' num2str(bloc) ' of ' num2str(nBloc) ': computing']); end
        % Project to taper space -> compute fft -> convert to power ->
        % average across tapers -> sqrt for consistency.
        J{bloc}=sqrt(mean(abs(reshape(fft(reshape(data{bloc}.*permute(tapers,[4 1 5 2 3]),[NT Cdist(bloc)*K]),NFFT)/Fs,[NFFT Cdist(bloc) K])).^2,3));
        toc

        % Delete data to save space
        data{bloc} = {};
        if verbose; disp(['bloc ' num2str(bloc) ' of ' num2str(nBloc) ': done']); end
        toc
    end
    % Catenate compute blocs
    J = cat(2,J{:});
end


%% Adjust spectrum and keep data with Fpass
df=Fs/NFFT;
f=0:df:(Fs-df);
f = f-Fs/2;
f = -f;
f = fftshift(f);
fInd = find(f>=Fpass(1) & f<=Fpass(2));
[~,b] = sort(f(fInd));
fInd = fInd(b);

sz = [length(fInd) K C/NR NR];
if avTapers; sz(2) = 1; end
J = reshape(permute(J(fInd,:,:),[1 3 2]),sz);
f = f(fInd);