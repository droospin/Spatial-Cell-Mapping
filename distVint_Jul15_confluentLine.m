
%%

% =======================
% Extract raw data
% =======================
distVint = track4out_Jan29(1).distVintIntProt1;

% =======================
% Remove rows that are all zeros
% =======================
distVint_clean = distVint(~all(distVint == 0, 2), :);

% =======================
% Plot metadata
% =======================
titleStg = 'Normalized F-actin vs. Cav1 Intensity';
xAxisStg = 'F-actin';
yAxisStg = 'Cav1';

lineXlim = [0 1];
lineYlim = [0 1];

Xtick = 0:0.2:1;
Ytick = 0:0.2:1;

% =======================
% Smooth + plot cleaned data
% =======================
interpDvI = smooth_and_plot( ...
    distVint_clean, ...
    0.02, ...
    titleStg, ...
    xAxisStg, ...
    yAxisStg, ...
    lineXlim, ...
    lineYlim, ...
    Xtick, ...
    Ytick);

% =======================
% Store outputs
% =======================
outDistVint = struct();
outDistVint.distVint  = distVint_clean;
outDistVint.interpDvI = interpDvI;
%% 

for i = 1:numel(track4out_Jan29)

    refPatt = double(track4out_Jan29(i).refProt);
    refPatt(track4out_Jan29(i).pattMask == 0) = 0;
    lineUP(i).refPatt = refPatt / max(refPatt(:));

    intPatt = double(track4out_Jan29(i).intProt1);
    intPatt(track4out_Jan29(i).pattMask == 0) = 0;
    lineUP(i).intPatt = intPatt / max(intPatt(:));

end

pattMaskVec       = [];
refProtPattVec    = [];
intProt1PattVec   = [];

for i = 1:numel(track4out_Jan29)
    pattMaskVec       = [pattMaskVec;     track4out_Jan29(i).pattMask(:)];
    refProtPattVec    = [refProtPattVec;  lineUP(i).refPatt(:)];
    intProt1PattVec   = [intProt1PattVec; lineUP(i).intPatt(:)];

end

lineUP(1).refProtPattVec  = refProtPattVec;
lineUP(1).intProt1PattVec = intProt1PattVec;

% RefProt vs. intProt1 (intensity only)
lineUP(1).refVSintPattVec(:, 1) = refProtPattVec;
lineUP(1).refVSintPattVec(:, 2) = intProt1PattVec;

%%

refVSint = track4out_Jan29(1).refVSintPattVec;
track(1).refVSint = refVSint(~all(refVSint == 0, 2), :);

titleStg       = 'Normalized F-actin vs. Cav1 Intensity';
xAxisStg       = 'F-actin';
yAxisStg       = 'Cav1';
lineXlim       = [0 1]; Xtick = 0:0.2:1;
lineYlim       = [0 1]; Ytick = 0:0.2:1;

[outRefVSint1] = smooth_and_plot(refVSint, 0.02, titleStg, xAxisStg, yAxisStg, lineXlim, lineYlim, Xtick, Ytick);

%%

% Initialize output vector
refProtCellsVec = [];

for i = 1:numel(track24out_Jan29)

    % Pull required fields
    refProtThresh = track24out_Jan29(i).refProt;
    normDistIm    = track24out_Jan29(i).intProt1normDist;  % same cell-defining mask

    % Normalize refProt per image (same philosophy as intProt1cells)
    refNorm = double(refProtThresh);
    maxVal  = double(max(refNorm(:)));

    if maxVal > 0
        refNorm = refNorm ./ maxVal;
    else
        continue;  % skip empty images safely
    end

    % Keep only traced cell pixels
    refNorm(normDistIm == 0) = 0;

    % Vectorize and remove zeros
    vals = refNorm(:);
    vals = vals(vals > 0);

    % Append
    refProtCellsVec = [refProtCellsVec; vals];

end

% Store once for downstream use
track(1).refProtCellsVec = refProtCellsVec;

distVint_clean(:, 1) = track24out_Jan29(1).distVintIntProt1(:, 1);
distVint_clean(:, 2) = refProtCellsVec(:, 1);

% =======================
% Plot metadata
% =======================
titleStg = 'Normalized Distance vs. F-actin Intensity';
xAxisStg = 'Distance';
yAxisStg = 'F-actin';

lineXlim = [0 1];
lineYlim = [0 1];

Xtick = 0:0.2:1;
Ytick = 0:0.2:1;

% =======================
% Smooth + plot cleaned data
% =======================

interpDvI = smooth_and_plot( ...
    distVint_clean, ...
    0.02, ...
    titleStg, ...
    xAxisStg, ...
    yAxisStg, ...
    lineXlim, ...
    lineYlim, ...
    Xtick, ...
    Ytick);


% =======================
% Store outputs
% =======================

track24_Jan29.distVintRef  = distVint_clean;
track24_Jan29.refInterpDvI = interpDvI;