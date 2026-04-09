function circleLocker = circle_3d(folder)
% circle_alignment:

%% Section 1: Get list of all TIF files in the folder, including subfolders

% Get list of all TIF files in the folder, including subfolders
imFiles = fullfile(folder, '**', '*.tif');
full    = dir(imFiles);

%% Section 2: Initialize Data Structure to store image results

%{
% circleLocker(length(full)) = struct('filename', [], ...
                                    'dna', [], 'refProt', [], 'intProt', [], ...
                                    'intProttranslate', [], 'intProtseg', [], ...
                                    'ogCent', [], 'translateCent', [], ...
                                    'refProtThresh', [], 'dnaThresh2x', []);
%}

numColors    = 256;                                                              % Define number of color levels
origColormap = hot(numColors);                                                   % Get the default 'hot' colormap
cutoff       = 0.9;                                                              % Define a cutoff threshold (e.g., 80% of the colormap)
newMaxIdx    = round(cutoff * numColors); 
customHot    = origColormap(1:newMaxIdx, :);                                     % Stretch and cut off colormap
customHot    = [customHot; repmat(customHot(end, :), numColors - newMaxIdx, 1)]; % Extend the last color to fill up to 1

imInterval = 1:512;
%% Section 3: Main loop; Process Each Image

%% -----------------------------
%% 0) List all files (cells) in 'full'
%% -----------------------------
% Assume “full” is already a struct array (from dir) containing all of your 
% TIFFs or MATs.  
numCells = length(full);

%% -----------------------------
%% 1) First pass: determine how many slices each cell has
%% -----------------------------
cellDepths = nan(numCells,1);
for i = 1:numCells
    source = fullfile(full(i).folder, full(i).name);
    data   = three_channel_stacks(source);
    % assume refProt stack gives the canonical z‐dimension
    cellDepths(i) = size(data.stackRefProt, 3);
end

% Decide on a common number of “interpolated z‐points” (Ntarget).
% You could choose something like the maximum depth across all cells,
% or a fixed value (e.g. 50). Here we’ll use max depth to preserve detail:
Ntarget = max(cellDepths);      % ← NORMALIZATION STEP: common grid length

% Preallocate big arrays: each cell will be resampled to [Ntarget × Nbins].
% We already know your radial‐distance binning uses 100 bins (0.01:0.01:1.00).
fixedDistances = (0.01:0.01:1.00)';   % 100×1 column vector  
numBins = numel(fixedDistances);

% These 3D arrays will store each cell’s interpolated “slice × distance” data:
allRefMean      = nan(numCells, Ntarget, numBins);   % one layer per cell
allIntMean      = nan(numCells, Ntarget, numBins);
allIntNoDNAmean = nan(numCells, Ntarget, numBins);
%% -----------------------------
%% 2) Second pass: process each cell and interpolate onto the common z‐grid
%% -----------------------------
for i = 1:numCells
    %% ---- read & segment exactly as before ----
    source = fullfile(full(i).folder, full(i).name);

    disp(full(i).name) % display file name
    
    data               = three_channel_stacks(source);
    dna                = data.stackDNA;    
    refProt            = data.stackRefProt;
    maxRefProt         = data.refProt;
    intProt            = data.stackIntProt;
    [xDim, yDim, zDim] = size(dna);

    dnaSeg = zeros(size(dna));
    for x = 1:xDim
        for y = 1:yDim
            for z = 1:zDim
                if dna(x, y, z) > 5500
                    dnaSeg(x, y, z) = dna(x, y, z);
                end
            end
        end
    end

    dnaBinSeg      = dnaSeg > 0;
    dnaRegionprops = regionprops(dnaBinSeg, 'Area');
    validRegions   = [dnaRegionprops.Area] >= 2000;
    dnaMask        = ismember(labelmatrix(bwconncomp(dnaBinSeg)), find(validRegions));

    dnaThresh2x    = zeros(size(dna));

    for x = 1:xDim
        for y = 1:yDim
            for z = 1:zDim
                if dnaMask(x, y, z) > 0
                    dnaThresh2x(x, y, z) = dnaSeg(x, y, z);
                end
            end
        end
    end

    %— [RefProt segmentation / thresholding] —
    refProtSeg = zeros(size(refProt));
    for x = 1:xDim
        for y = 1:yDim
            for z = 1:zDim
                if refProt(x, y, z) > 4000
                    refProtSeg(x, y, z) = refProt(x, y, z);
                end
            end
        end
    end

    %— [IntProt segmentation / thresholding] —
    intProtSeg      = zeros(size(intProt));

    for x = 1:xDim
        for y = 1:yDim
            for z = 1:zDim
                if intProt(x, y, z) > 4000
                   intProtSeg(x, y, z) = intProt(x, y, z);
                end
            end
        end
    end

    intProtSegNoDNA = intProtSeg;

    for x = 1:xDim
        for y = 1:yDim
            for z = 1:zDim
                if dnaThresh2x(x, y, z) > 0
                   intProtSegNoDNA(x, y, z) = 0;
                end
            end
        end
    end

    %% Section 4: Compute Image Contrast & Adjust Segmentation Sensitivity
    
    contrast            = std(double(refProt(:))) / mean(double(refProt(:)));         % Compute image contrast metric
    adaptiveSensitivity = max(0.3, 0.5 - 0.2 * contrast);                             % Adjust sensitivity dynamically (tweak scaling factor as needed)
    maxRefProtThresh    = segment_singlecell_circle(maxRefProt, adaptiveSensitivity); % Call the segmentation function with desired sensitivity

    %% Section 5: Create Binary Mask & Extract Centroid
    
    refProtMask = maxRefProtThresh > 0; % Convert to binary mask

    refProtRegionprops = regionprops(refProtMask, 'Centroid', 'Area'); % Find the centroid of the current image

    if isempty(refProtRegionprops)                                         % If no regions were found,
       warning('No refProt regions found in %s. Skipping.', full(i).name); % then display message
       continue;                                                           % and skip this image.
    end
    
    % Find the centroid of the largest region just in case there are multiple.
    [~, largestIdx] = max([refProtRegionprops.Area]); 
    cellCent        = refProtRegionprops(largestIdx).Centroid;
        
    %% Section 6: Define Reference Centroid

    % If this is the first image, set the reference centroid
    if i == 1
        imCent = [256, 256];
    end

    translateVec          = imCent - cellCent;
    refProttranslate      = imtranslate(refProtSeg, translateVec);
    intProttranslate      = imtranslate(intProtSeg, translateVec);
    intProtNoDNAtranslate = imtranslate(intProtSegNoDNA, translateVec);

    %— Clean up RefProt slice‐by‐slice (area ≥100) as before —
    refProtMask   = refProttranslate > 0;
    refProtThresh = zeros(size(refProttranslate));

    for z = 1:size(refProttranslate, 3)
        sliceMask          = refProtMask(:,:,z);
        cc                 = bwconncomp(sliceMask);
        stats              = regionprops(cc, 'Area');
        validIdx           = find([stats.Area] >= 100);
        cleanMask          = ismember(labelmatrix(cc), validIdx);
        refProtMask(:,:,z) = cleanMask;
    
        slice              = double(refProttranslate(:,:,z));    % <-- convert to double
        mask               = double(cleanMask);                  % <-- also convert mask
        refProtThresh(:,:,z) = slice .* mask;                    % <-- safe to multiply
    end


    %— Compute Euc Dist from center (256,256) for all positive pixels in each slice —
    REFeucDistCent      = zeros(size(refProt));
    INTeucDistCent      = zeros(size(intProt));
    INTNoDNAeucDistCent = zeros(size(intProt));

    for x = imInterval
        for y = imInterval
            for z = 1:size(refProt, 3)
                if refProtThresh(x, y, z) > 0
                   REFeucDistCent(x, y, z) = sqrt((256 - x)^2 + (256 - y)^2);
                end
                if intProttranslate(x, y, z) > 0
                   INTeucDistCent(x, y, z) = sqrt((256 - x)^2 + (256 - y)^2);
                end
                if intProtNoDNAtranslate(x, y, z) > 0
                   INTNoDNAeucDistCent(x, y, z) = sqrt((256 - x)^2 + (256 - y)^2);
                end
            end
        end
    end

    %% -----------------------------
    %% 2a) Normalize radial distances and intensities (unchanged)
    %% -----------------------------
    maxRefEucDist       = max(REFeucDistCent(:));
    normIntEucDist      = round(INTeucDistCent / maxRefEucDist, 2);
    normRefEucDist      = round(REFeucDistCent / maxRefEucDist, 2);
    normIntNoDNAeucDist = round(INTNoDNAeucDistCent / maxRefEucDist, 2);

    normRefInt      = refProtThresh / max(refProtThresh(:));
    normIntInt      = intProttranslate / max(intProttranslate(:));
    normIntNoDNAint = intProtNoDNAtranslate / max(intProtNoDNAtranslate(:));

    %% -----------------------------
    %% 2b) Build “slice × pixel‐vector” for distance and intensity
    %% -----------------------------
    Zi          = size(refProt,3);  % number of raw slices for this cell
    M           = size(refProt,1)*size(refProt,2);

    refDistVec      = zeros(Zi, M);
    intDistVec      = zeros(Zi, M);
    intNoDNAdistVec = zeros(Zi, M);

    refIntVec      = zeros(Zi, M);
    intIntVec      = zeros(Zi, M);
    intNoDNAintVec = zeros(Zi, M);

    for z = 1:Zi
        refDistVec(z, :)      = reshape(normRefEucDist(:, :, z), 1, []);
        intDistVec(z, :)      = reshape(normIntEucDist(:, :, z), 1, []);
        intNoDNAdistVec(z, :) = reshape(normIntNoDNAeucDist(:,:,z), 1, []);

        refIntVec(z, :)  = reshape(normRefInt(:, :, z),  1, []);
        intIntVec(z, :)  = reshape(normIntInt(:, :, z),  1, []);
        intNoDNAintVec(z, :) = reshape(normIntNoDNAint(:, :, z), 1, []);
    end

    %% -----------------------------
    %% 2c) For each slice, bin the intensities by fixedDistances
    %% -----------------------------
    %% -----------------------------
%% 2c) For each slice, bin the intensities by fixedDistances
%% -----------------------------
% Define bin edges once, outside the loop
binEdges = [fixedDistances - 0.005; 1.005];  % Creates 100 bins spanning 0.01 to 1.00

refMeanInt_i      = nan(Zi, numBins);
intMeanInt_i      = nan(Zi, numBins);
intNoDNAmeanInt_i = nan(Zi, numBins);

for z = 1:Zi
    d_ref = refDistVec(z, :)';
    i_ref = refIntVec(z, :)';

    d_int = intDistVec(z, :)';
    i_int = intIntVec(z, :)';

    d_intNoDNA = intNoDNAdistVec(z, :)';
    i_intNoDNA = intNoDNAintVec(z, :)';

    % Bin distances robustly using discretize (avoids floating-point errors)
    idxRef_clean      = discretize(d_ref, binEdges);
    idxInt_clean      = discretize(d_int, binEdges);
    idxIntNoDNA_clean = discretize(d_intNoDNA, binEdges);

    % Keep only valid indices (non-NaN)
    validRef      = ~isnan(idxRef_clean);
    validInt      = ~isnan(idxInt_clean);
    validIntNoDNA = ~isnan(idxIntNoDNA_clean);

    i_ref_clean      = i_ref(validRef);
    i_int_clean      = i_int(validInt);
    i_intNoDNA_clean = i_intNoDNA(validIntNoDNA);

    idxRef_clean      = idxRef_clean(validRef);
    idxInt_clean      = idxInt_clean(validInt);
    idxIntNoDNA_clean = idxIntNoDNA_clean(validIntNoDNA);

    % Accumulate mean intensity per bin
    tmpRef      = accumarray(idxRef_clean, i_ref_clean, [numBins,1], @mean, NaN);
    tmpInt      = accumarray(idxInt_clean, i_int_clean, [numBins,1], @mean, NaN);
    tmpIntNoDNA = accumarray(idxIntNoDNA_clean, i_intNoDNA_clean, [numBins, 1], @mean, NaN);

    refMeanInt_i(z, :)      = tmpRef';
    intMeanInt_i(z, :)      = tmpInt';
    intNoDNAmeanInt_i(z, :) = tmpIntNoDNA';
end

    %% -----------------------------
    %% 2d) INTERPOLATE each cell’s (Zi × numBins) onto (Ntarget × numBins)
    %% (Option 1: normalize z‐axis).
    %% -----------------------------
    % Build normalized z‐coordinate for this cell:
    % If Zi == 1, we’ll just replicate that row across all targets.
    if Zi == 1
        t_orig = 0;
    else
        t_orig = linspace(0, 1, Zi);
    end

    t_target = linspace(0, 1, Ntarget);

    % Preallocate
    refInterp      = nan(Ntarget, numBins);
    intInterp      = nan(Ntarget, numBins);
    intNoDNAInterp = nan(Ntarget, numBins);

    for b = 1:numBins
        % Take column b of refMeanInt_i: that is a Zi×1 vector
        v_ref      = refMeanInt_i(:, b);
        v_int      = intMeanInt_i(:, b);
        v_intNoDNA = intNoDNAmeanInt_i(:, b);

        % If Zi == 1, interp1 will fail. Just tile the single row:
        if Zi == 1
            refInterp(:, b)      = repmat(v_ref, Ntarget, 1);
            intInterp(:, b)      = repmat(v_int, Ntarget, 1);
            intNoDNAInterp(:, b) = repmat(v_intNoDNA, Ntarget, 1);
        else
            % Use linear interpolation; “extrap” set to NaN by default if out of range
            refInterp(:, b)      = interp1(t_orig, v_ref, t_target, 'linear', NaN);
            intInterp(:, b)      = interp1(t_orig, v_int, t_target, 'linear', NaN);
            intNoDNAInterp(:, b) = interp1(t_orig, v_intNoDNA, t_target, 'linear', NaN);
        end
    end

    % Store into the big 3D accumulator (cell i’s layer)
    allRefMean(i,:,:)        = refInterp;  
    allIntMean(i,:,:)        = intInterp;
    allIntNoDNAmean(i, :, :) = intNoDNAInterp;

    %% Section 7: Store Original & Processed Data

    circleLocker(i).filename = full(i).name;

    circleLocker(i).refProt = refProt;
    circleLocker(i).intProt = intProt;
    circleLocker(i).dna     = dna;

    circleLocker(i).refProttranslate   = refProttranslate;
    circleLocker(i).intProtThresh      = intProttranslate;
    circleLocker(i).intProtNoDNAthresh = intProtNoDNAtranslate;

    circleLocker(i).refProtThresh    = refProtThresh;
    circleLocker(i).dnaThresh        = dnaThresh2x;

    circleLocker(i).imCent      = imCent;
    circleLocker(i).initialCent = cellCent;

    circleLocker(i).rawRefSideProfile      = refMeanInt_i;      % size Zi×numBins
    circleLocker(i).rawIntSideProfile      = intMeanInt_i;      % size Zi×numBins
    circleLocker(i).rawIntNoDNASideProfile = intNoDNAmeanInt_i;

    circleLocker(i).allRefMean      = allRefMean;
    circleLocker(i).allIntMean      = allIntMean;
    circleLocker(i).allIntNoDNAmean = allIntNoDNAmean;

    circleLocker(i).depth = cellDepths(i);

    circleLocker(i).fixedDistances = fixedDistances;
    circleLocker(i).Ntarget        = Ntarget;
    circleLocker(i).t_target       = linspace(0,1,Ntarget);

    circleLocker(i).thresholds.dna        = 5500;
    circleLocker(i).thresholds.refProt    = 3500;
    circleLocker(i).thresholds.minAreaDNA = 2000;
    circleLocker(i).thresholds.minAreaRef = 100;
end

%% -----------------------------
%% 3) After the loop: average across cells (ignore NaNs)
%% -----------------------------
% allRefMean is [numCells × Ntarget × numBins]
% We want the mean over the first dimension:
avgRefMean      = squeeze(nanmean(allRefMean, 1));           % size: [Ntarget × numBins]
avgIntMean      = squeeze(nanmean(allIntMean, 1));           % size: [Ntarget × numBins]
avgIntNoDNAmean = squeeze(nanmean(allIntNoDNAmean, 1));

avgRefMean      = avgRefMean / max(avgRefMean(:));
avgIntMean      = avgIntMean / max(avgIntMean(:));
avgIntNoDNAmean = avgIntNoDNAmean / max(avgIntNoDNAmean(:));

%% -----------------------------
%% 4) Plot the final “mean across cells” images
%% -----------------------------

figure;
imagesc(fixedDistances, linspace(1, Ntarget, Ntarget), avgRefMean);
set(gca, 'YDir','normal');
xlabel('Normalized Distance (bin)');
ylabel('Normalized Slice Index');
%title('⟨Reference Channel⟩ – Averaged Over All Cells');
colormap(hot);
%colorbar;
axis off;

figure;
imagesc(fixedDistances, linspace(1, Ntarget, Ntarget), avgIntMean);
set(gca, 'YDir','normal');
xlabel('Normalized Distance (bin)');
ylabel('Normalized Slice Index');
%title('⟨Interest Channel⟩ – Averaged Over All Cells');
colormap(hot);
%colorbar;
axis off;

figure;
imagesc(fixedDistances, linspace(1, Ntarget, Ntarget), avgIntNoDNAmean);
set(gca, 'YDir','normal');
xlabel('Normalized Distance (bin)');
ylabel('Normalized Slice Index');
%title('⟨Interest Channel w/o Nuclear Signal⟩ – Averaged Over All Cells');
colormap(hot);
%colorbar;
axis off;

%% Save averaged info

circleLocker(1).avgRefMean      = avgRefMean;
circleLocker(1).avgIntMean      = avgIntMean;
circleLocker(1).avgIntNoDNAmean = avgIntNoDNAmean;

end % end main function

%%

function imThresh = segment_singlecell_circle(imX, sensitivity)
% segment_singlecell_circle: Segments a single cell using adaptive thresholding
% Generated with help from the Image Segmenter app on 17-Feb-2025
%----------------------------------------------------


imX    = imadjust(imX); % Adjust data to span full intensity range

BWim   = imbinarize(im2gray(imX), 'adaptive', 'Sensitivity', sensitivity, 'ForegroundPolarity', 'dark'); % Apply adaptive thresholding for segmentation

BWfill = imfill(BWim, 'holes'); % Fill Holes in Binary Mask

% Remove small regions below a threshold size
BWregionprops = regionprops(BWfill, 'Area');
sizeThresh1   = 100000;
for i = 1:numel(BWregionprops)                                 
    if BWregionprops(i).Area < sizeThresh1
       BWregionprops(i).Area = 0;                               
    end
end
BWmask1 = ismember(labelmatrix(bwconncomp(BWfill)), find([BWregionprops.Area] > 0));

BWerode = imerode(BWmask1, strel('disk', 1)); % Erosion to Refine Mask

% Remove Small Objects
sizeThresh2 = 100;
BWregionprops = regionprops(BWerode, 'Area');
for i = 1:numel(BWregionprops)                                 
    if BWregionprops(i).Area < sizeThresh2
       BWregionprops(i).Area = 0;                               
    end
end
BWmask2  = ismember(labelmatrix(bwconncomp(BWerode)), find([BWregionprops.Area] > 0));

BWdilate = imdilate(BWmask2, strel('disk', 1)); % Dilation to restore eroded perimeter

% Final Size Thresholding
sizeThresh3 = 8000;
BWregionprops = regionprops(BWdilate, 'Area');
for i = 1:numel(BWregionprops)                                 
    if BWregionprops(i).Area < sizeThresh3
       BWregionprops(i).Area = 0;                               
    end
end
BWmask3 = ismember(labelmatrix(bwconncomp(BWdilate)), find([BWregionprops.Area] > 0));

% Create Masked Image
imThresh = imX;
imThresh(~BWmask3) = 0;

end