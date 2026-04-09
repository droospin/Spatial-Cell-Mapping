function outLocker = backSubtract_Renorm(locker, backThresh)

%% Presets
if nargin < 2
    backThresh = 5000;
end

outLocker = locker; % Set up variable

excludeFields = {'intProt1normDist'};

fn = fieldnames(locker);

for j = 1:numel(locker)

    maxInt        = max(locker(j).intProt1(:));
    backThreshDec = double(backThresh ./ maxInt);

for i = 1:numel(fn)

    fname      = fn{i};
    intensity  = locker(j).(fname);

    % Skip non-numeric or empty
    if ~isnumeric(intensity) || isempty(intensity)
       outLocker(j).(fname) = intensity;
       continue
    end

    if maxInt == 0
       outLocker(j).(fname) = intensity;
       continue
    end

    % Skip excluded categories
    if any(strcmp(fname, excludeFields))
       outLocker(j).(fname) = locker(j).(fname);
       continue
    end

    %% Non-normalized intensity variables
    if strcmp(fname, 'intProt1')
       
    intensity                = intensity - backThresh;
    intensity(intensity < 0) = 0;

    outLocker(j).(fname) = intensity;

    %% VS variables: Nx2, modify column 2 
    elseif contains(fname, 'VS')

    if size(intensity, 2) ~= 2
       warning('%s marked VS but not Nx2 — skipped', fname);
       outLocker(j).(fname) = intensity;
       continue
    end
    
    col2 = intensity(:, 2);
    
    if max(col2) == 0
       outLocker(j).(fname) = intensity;
       continue
    end
    
    disp(min(col2(:)));
    col2           = col2 - backThreshDec;
    col2(col2 < 0) = 0;
    col2           = col2 ./ max(col2);
    disp(min(col2(:)));

    intensity(:, 2)      = col2;
    outLocker(j).(fname) = intensity;

    %% Normalized intensity variables
    elseif contains(fname, 'intProt1') && ~contains(fname, 'VS') && ~strcmp(fname, 'intProt1')

    if max(intensity(:)) == 0
       outLocker(j).(fname) = intensity;
       continue
    end
    
    disp(min(intensity(:)));
    intensity                = intensity - backThreshDec;
    intensity(intensity < 0) = 0;

    maxVal = max(intensity(:));
    if maxVal > 0
       intensity = intensity ./ maxVal;
    end
    disp(min(intensity(:)));

    outLocker(j).(fname) = intensity;    

    %% for variables not considered but need to be kept (ex: Ref, etc)    
    else

    outLocker(j).(fname) = locker(j).(fname);    
    
    end % end if else 

end % end for loop
end % numel locker
end % end function

%{
'intProt1', 'intProt1cells', 'intProt1pattEdge',
'intProt1junc', 'intProt1edge', 'intProt1NOjunc',
'intProt1PattVec', 'intProt1cellsVec', 'intProt1pattEdgeVec',
'intProt1juncVec', 'intProt1edgeVec', 'intProt1NOjuncVec',
'refVSintPattVec', 'distVintIntProt1'
%}