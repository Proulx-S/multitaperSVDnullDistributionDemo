function perms = replacePerm(perms,perm,whichEnd,validate)
% perms is meant to store randomization values, where the last dimension is
% contains n realizations of the randomization process.
% perm should be a new realisation of the randomization process.
% If a value of perm is larger than the smallest value of perms, the latter
% is replaced with the former.
% This is meant to save memory when empirically estimating a distribution,
% by storing only the distribution's extrem values.
if ~exist('validate','var'); validate = false; end
if validate; permsBefore = perms; end



permDim = length(size(perms));
switch whichEnd
    case 'top'
        [minInMx,index] = min(abs(perms),[],permDim);  %gives global index and a new matrix corresponding to those values
        idxToReplace = find(abs(perm)-minInMx > 0);  %global index
    case 'bot'
        [maxInMx,index] = max(abs(perms),[],permDim);  %gives global index and a new matrix corresponding to those values
        idxToReplace = find(maxInMx-abs(perm) > 0);  %global index
end
indStr = strjoin(cellstr(num2str((1:permDim-1)','ind%d'))',',');
eval(['[' indStr '] = ind2sub(size(perms),idxToReplace);']); % [xx,yy,zz] = ind2sub(size(perms),idxToReplace);
eval(['perms(sub2ind(size(perms),' indStr ',index(idxToReplace))) = perm(idxToReplace);']); % perms(sub2ind(size(perms),xx,yy,zz,index(idxToReplace))) = perm(idxToReplace);  %do the replacement







if validate
    %%% double-check the replacement does what it is supposed to
    nPerms = size(perms,permDim);

    N = nnz(permsBefore~=perms);
    % disp([num2str(N) ' values replaced'])

    if N
        % never replaces more than 1 value
        tmp1 = sum(permsBefore~=perms,4)>1;
        tmp1 = tmp1(:);
        any(tmp1); % should be 0
        if any(tmp1)
            warning('something wrong');
            keyboard
        end

        % replaced values are always replaced by a larger value
        switch whichEnd
            case 'top'
                % larger value
                tmp2 = abs(perms(permsBefore~=perms)) > abs(permsBefore(permsBefore~=perms));
            case 'bot'
                % smaller value
                tmp2 = abs(perms(permsBefore~=perms)) < abs(permsBefore(permsBefore~=perms));
        end
        all(tmp2); % should be 1
        if ~all(tmp2)
            warning('something wrong');
            keyboard
        end


        % replaced values are always
        switch whichEnd
            case 'top'
                % the smallest ones
                uTop5BeforeMin = repmat( min(abs(permsBefore),[],4) , [1 1 1 nPerms] );
                tmp3 = abs(permsBefore(permsBefore~=perms)) > abs(uTop5BeforeMin(permsBefore~=perms));
            case 'bot'
                % the largest ones
                uBot5BeforeMax = repmat( max(abs(permsBefore),[],4) , [1 1 1 nPerms] );
                tmp3 = abs(permsBefore(permsBefore~=perms)) < abs(uBot5BeforeMax(permsBefore~=perms));
        end
        any(tmp3); % should be 0
        if any(tmp3)
            warning('something wrong');
            keyboard
        end

        % disp([any(tmp1) ~all(tmp2) any(tmp3)]) % should all be 0
        % if any([any(tmp1) ~all(tmp2) any(tmp3)]); error('something wrong'); end
    end
end