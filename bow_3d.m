function [bowLocker] = bow_3d(folder)

%% Section 1: Get list of all TIF files in the folder, including subfolders

% Get list of all TIF files in the folder, including subfolders
imFiles = fullfile(folder, '**', '*.tif');
full    = dir(imFiles);

%% Section 2: Initialize Data Structure to store image results

numColors    = 256;                                                              % Define number of color levels
origColormap = hot(numColors);                                                   % Get the default 'hot' colormap
cutoff       = 0.9;                                                              % Define a cutoff threshold (e.g., 80% of the colormap)
newMaxIdx    = round(cutoff * numColors); 
customHot    = origColormap(1:newMaxIdx, :);                                     % Stretch and cut off colormap
customHot    = [customHot; repmat(customHot(end, :), numColors - newMaxIdx, 1)]; % Extend the last color to fill up to 1

imInterval     = 1:512;
numBins        = 100;                   % or whatever number of bins you want
fixedDistances = linspace(-1, 1, numBins);
Ntarget        = 50;                    % or your desired z-slices

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
fixedDistances = (-1.00: 0.01: 1.00)';   % 100×1 column vector  
numBins        = numel(fixedDistances);

% These 3D arrays will store each cell’s interpolated “slice × distance” data:
allRefMean      = nan(numCells, Ntarget, numBins);   % one layer per cell
allIntMean      = nan(numCells, Ntarget, numBins);

%% -----------------------------
%% 2) Second pass: process each cell and interpolate onto the common z‐grid
%% -----------------------------
for i = 1:numCells

    % If this is the first image, set the reference centroid
    if i == 1
        imCent = [256, 256];
    end
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

    dnaThresh    = zeros(size(dna));

    for x = 1:xDim
        for y = 1:yDim
            for z = 1:zDim
                if dnaMask(x, y, z) > 0
                    dnaThresh(x, y, z) = dnaSeg(x, y, z);
                end
            end
        end
    end

    %%

    %% Section 10: Compute Image Contrast & Adjust Segmentation Sensitivity
    
    contrast            = std(double(maxRefProt(:))) / mean(double(maxRefProt(:))); % Compute image contrast metric
    adaptiveSensitivity = max(0.1, 0.3 - 0.2 * contrast);                           % Adjust sensitivity dynamically (tweak scaling factor as needed)
    maxRefProtThresh    = segment_singlecell_bow(maxRefProt, adaptiveSensitivity);  % Call the segmentation function with desired sensitivity
    
    imshow(maxRefProtThresh > 0);
    hold on;
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

    plot(cellCent(1), cellCent(2), 'rx', 'MarkerSize', 10, 'LineWidth', 2); pause (0.5);
    hold off;

    %%

    %— [RefProt segmentation / thresholding] —
    refProtSeg = zeros(size(refProt));
    for x = 1:xDim
        for y = 1:yDim
            for z = 1:zDim
                if maxRefProtThresh(x, y) > 0
                   refProtSeg(x, y, z) = refProt(x, y, z);
                end
            end
        end
    end

    %— [RefProt segmentation / thresholding] —
    for x = 1:xDim
        for y = 1:yDim
            for z = 1:zDim
                if refProtSeg(x, y, z) < 3500
                   refProtSeg(x, y, z) = 0;
                end
            end
        end
    end

    %%
    
    %— [IntProt segmentation / thresholding] —
    intProtSeg = zeros(size(intProt));

    for x = 1:xDim
        for y = 1:yDim
            for z = 1:zDim
                if maxRefProtThresh(x, y) > 0
                   intProtSeg(x, y, z) = intProt(x, y, z);
                end
            end
        end
    end

    for x = 1:xDim
        for y = 1:yDim
            for z = 1:zDim
                if intProtSeg(x, y, z) < 3500
                   intProtSeg(x, y, z) = 0;
                end
            end
        end
    end

    imshow(intProtSeg(:, :, 7), []); pause (0.5);

    %% Section 6: Define Reference Centroid
    
    translateVec     = imCent - cellCent;
    refProttranslate = zeros(size(refProt));
    intProttranslate = zeros(size(refProt));

    for z = 1:size(refProt, 3)
        refProttranslate(:, :, z) = imtranslate(refProtSeg(:, :, z), translateVec);
        intProttranslate(:, :, z) = imtranslate(intProtSeg(:, :, z), translateVec);
    end

    imshow(refProttranslate(:, :, 7), []); pause(0.5); 

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
    
        slice                = double(refProttranslate(:,:,z));  % <-- convert to double
        mask                 = double(cleanMask);                % <-- also convert mask
        refProtThresh(:,:,z) = slice .* mask;                    % <-- safe to multiply
    end


    %— Compute Euc Dist from center (256,256) for all positive pixels in each slice —
    REFeucDistCent = nan(xDim, yDim, zDim);
    
    for x = 1:xDim
        for y = 1:yDim
            for z = 1:zDim
                if refProtThresh(x, y, z) > 0
                   if x < 256
                      dist = sqrt((x - imCent(1))^2 + (y - imCent(2))^2);
                   else
                      dist = - sqrt((x - imCent(1))^2 + (y - imCent(2))^2);
                   end

                   REFeucDistCent(x, y, z) = dist;

                end
            end
        end
    end

    imshow(REFeucDistCent(:, :, 7), []); pause(0.5);

    %{
    %% For cropping -- Section Optional!
    %— Compute Euc Dist from center (256,256) for all positive pixels in each slice —
    REFDistY = nan(xDim, yDim, zDim);
    
    for x = 1:xDim
        for y = 1:yDim
            for z = 1:zDim
                if refProttranslate(x, y, z) > 0

                      dist              = abs(y - imCent(2));
                      REFDistY(x, y, z) = dist;

                end
            end
        end
    end

    % normalize Y variation
    maxY = max(REFDistY(~isnan(REFDistY)));

    REFDistY(~isnan(REFDistY)) = REFDistY(~isnan(REFDistY)) / maxY;

    imshow(REFDistY(:,:, 7), []); pause(0.5);

    for x = 1:xDim
        for y = 1:yDim
            for z = 1:zDim
                if REFDistY(x, y, z) > 0.2 || refProtThresh(x, y, z) == 0

                   REFeucDistCent(x, y, z) = nan;

                end
            end
        end
    end

    imshow(REFeucDistCent(:, :, 7), []); pause(0.5);
    %}

    %% -----------------------------
    %% 2a) Normalize radial distances and intensities (unchanged)
    %% -----------------------------
    maxRefEucDist  = max(REFeucDistCent(~isnan(REFeucDistCent)));

    % initialize
    normRefEucDist = nan(size(REFeucDistCent));

    normRefInt = zeros(size(refProtThresh));
    normIntInt = zeros(size(intProttranslate));

    % annotate valid pixels
    validMask = refProtThresh > 0  & ~isnan(REFeucDistCent);


    normRefEucDist(validMask) = REFeucDistCent(validMask) / maxRefEucDist;

    normRefInt(validMask) = refProtThresh(validMask)    / max(refProtThresh(:));
    normIntInt(validMask) = intProttranslate(validMask) / max(intProttranslate(:));

    sum(normRefEucDist(:) == 0)

    %% -----------------------------
    %% 2b) Build “slice × pixel‐vector” for distance and intensity
    %% -----------------------------
    Zi = size(refProt, 3);  % number of raw slices for this cell
    M  = 512 * 512;

    refDistVec = zeros(Zi, M);

    refIntVec = zeros(Zi, M);
    intIntVec = zeros(Zi, M);

    for z = 1:Zi
        refDistVec(z, :) = reshape(normRefEucDist(:, :, z), 1, []);

        refIntVec(z, :) = reshape(normRefInt(:, :, z),  1, []);
        intIntVec(z, :) = reshape(normIntInt(:, :, z),  1, []);
    end

    %% -----------------------------
    %% 2c) For each slice, bin the intensities by fixedDistances
    %% -----------------------------
    % Define bin edges once, outside the loop
    binEdges = [fixedDistances - 0.005; 1.005];  % Creates 100 bins spanning 0.01 to 1.00
    
    refMeanInt_i = nan(Zi, numBins);
    intMeanInt_i = nan(Zi, numBins);
    
    for z = 1:Zi
        d_ref = refDistVec(z, :)';
        i_ref = refIntVec(z, :)';
    
        d_int = refDistVec(z, :)';
        i_int = intIntVec(z, :)';
    
        % Bin distances robustly using discretize (avoids floating-point errors)
        idxRef_clean = discretize(d_ref, binEdges);
        idxInt_clean = discretize(d_int, binEdges);
    
        % Keep only valid indices (non-NaN)
        validRef = ~isnan(idxRef_clean);
        validInt = ~isnan(idxInt_clean);
    
        i_ref_clean = i_ref(validRef);
        i_int_clean = i_int(validInt);
    
        idxRef_clean = idxRef_clean(validRef);
        idxInt_clean = idxInt_clean(validInt);
    
        % Accumulate mean intensity per bin
        tmpRef = accumarray(idxRef_clean, i_ref_clean, [numBins,1], @mean, NaN);
        tmpInt = accumarray(idxInt_clean, i_int_clean, [numBins,1], @mean, NaN);
    
        refMeanInt_i(z, :) = tmpRef';
        intMeanInt_i(z, :) = tmpInt';
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
    refInterp = nan(Ntarget, numBins);
    intInterp = nan(Ntarget, numBins);

    for b = 1:numBins
        % Take column b of refMeanInt_i: that is a Zi×1 vector
        v_ref = refMeanInt_i(:, b);
        v_int = intMeanInt_i(:, b);

        % If Zi == 1, interp1 will fail. Just tile the single row:
        if Zi == 1
            refInterp(:, b) = repmat(v_ref, Ntarget, 1);
            intInterp(:, b) = repmat(v_int, Ntarget, 1);
        else
            % Use linear interpolation; “extrap” set to NaN by default if out of range
            refInterp(:, b) = interp1(t_orig, v_ref, t_target, 'linear', NaN);
            intInterp(:, b) = interp1(t_orig, v_int, t_target, 'linear', NaN);
        end
    end

    % Store into the big 3D accumulator (cell i’s layer)
    allRefMean(i,:,:) = refInterp;  
    allIntMean(i,:,:) = intInterp;

        %% Section 7: Store Original & Processed Data

    bowLocker(i).filename = full(i).name;

    bowLocker(i).refProt = refProt;
    bowLocker(i).intProt = intProt;
    bowLocker(i).dna     = dna;

    bowLocker(i).intProtSeg = intProtSeg;
    bowLocker(i).refProtSeg = refProtSeg;

    bowLocker(i).refProttranslate = refProttranslate;
    bowLocker(i).intProttranslate = intProttranslate;

    bowLocker(i).refProtThresh = refProtThresh;
    bowLocker(i).dnaThresh     = dnaThresh;

    bowLocker(i).imCent      = imCent;
    bowLocker(i).initialCent = cellCent;

    bowLocker(i).REFeucDistCent = REFeucDistCent;

    bowLocker(i).rawRefSideProfile = refMeanInt_i; % size Zi×numBins
    bowLocker(i).rawIntSideProfile = intMeanInt_i; % size Zi×numBins

    bowLocker(i).allRefMean = allRefMean;
    bowLocker(i).allIntMean = allIntMean;

    bowLocker(i).depth = cellDepths(i);

    bowLocker(i).fixedDistances = fixedDistances;
    bowLocker(i).Ntarget        = Ntarget;
    bowLocker(i).t_target       = linspace(0,1,Ntarget);

    bowLocker(i).thresholds.dna        = 5500;
    bowLocker(i).thresholds.refProt    = 3500;
    bowLocker(i).thresholds.minAreaDNA = 2000;
    bowLocker(i).thresholds.minAreaRef = 100;

end

    %% -----------------------------
    %% 3) After the loop: average across cells (ignore NaNs)
    %% -----------------------------
    % allRefMean is [numCells × Ntarget × numBins]
    % We want the mean over the first dimension:
    avgRefMean      = squeeze(nanmean(allRefMean, 1));           % size: [Ntarget × numBins]
    avgIntMean      = squeeze(nanmean(allIntMean, 1));           % size: [Ntarget × numBins]
    
    avgRefMean      = avgRefMean / max(avgRefMean(:));
    avgIntMean      = avgIntMean / max(avgIntMean(:));

    bowLocker(1).avgRefMean = avgRefMean;
    bowLocker(1).avgIntMean = avgIntMean;
    
    %% -----------------------------
    %% 4) Plot the final “mean across cells” images
    %% -----------------------------
    
    figure;
    imagesc(fixedDistances, linspace(1, Ntarget, Ntarget), avgRefMean);
    set(gca, 'YDir','normal');
    %xlabel('Normalized Distance (bin)');
    %ylabel('Normalized Slice Index');
    %title('⟨Reference Channel⟩ – Averaged Over All Cells');
    colormap(hot);
    %colorbar;
    axis off;

    figure;
    imagesc(fixedDistances, linspace(1, Ntarget, Ntarget), avgIntMean);
    set(gca, 'YDir','normal');
    %xlabel('Normalized Distance (bin)');
    %ylabel('Normalized Slice Index');
    %title('⟨Interest Channel⟩ – Averaged Over All Cells');
    colormap(hot);
    %colorbar;
    axis off;

    % end main loop
end

%% Section 13: Helper Function to rotate images based on orientation to x-axis in raw images

function varargout = rotate_and_crop(angle, padSize, varargin)
% rotate_and_crop: Rotates multiple images while preserving structure at edges.
% Usage:
%   [out1, out2, ...] = rotate_and_crop(angle, padSize, in1, in2, ...)
% 
% Inputs:
%   - angle: Rotation angle in degrees
%   - padSize: Padding size (in pixels)
%   - varargin: Input images to process
%
% Outputs:
%   - varargout: Rotated and cropped images
    
numImages = length(varargin);   % Number of images passed
varargout = cell(1, numImages); % Prepare output cell array
    
for n = 1:numImages
    
    padIm        = padarray(varargin{n}, [padSize, padSize], 'replicate', 'both'); % Pad image to prevent loss at edges
    rotateIm     = imrotate(padIm, angle, 'bilinear', 'crop');                     % Rotate the padded image
    ogCropIm     = rotateIm(padSize+1:end-padSize, padSize+1:end-padSize);         % Crop back to original size
    varargout{n} = ogCropIm;                                                       % Store the result

end
end

%% Section 14: Helper function to threshold cell mask channel

function imThresh = segment_singlecell_bow(imX, sensitivity)
% Segment image using auto-generated code from Image Segmenter app
% [BW,MASKEDIMAGE] = segmentImage(X) segments image X using auto-generated
% code from the Image Segmenter app. The final segmentation is returned in
% BW, and a masked image is returned in MASKEDIMAGE.
    
% Partly auto-generated by imageSegmenter app on 17-Feb-2025
%----------------------------------------------------

imX = imadjust(imX); % Adjust data to span data range
    
BWim = imbinarize(im2gray(imX), 'adaptive', 'Sensitivity', sensitivity, 'ForegroundPolarity', 'dark'); % Threshold image with adaptive threshold

BWregionprops = regionprops(BWim, 'Area', 'BoundingBox'); % Compute properties of connected components in the binary image BWfill

minSizeThresh = 4000;  % Define minimum
maxSizeThresh = 35000; % & maximum size thresholds for filtering regions
imHeight      = 512;

% Loop through each detected region
for i = 1:numel(BWregionprops)                                 
    if BWregionprops(i).Area < minSizeThresh || BWregionprops(i).Area > maxSizeThresh % If the region's area is outside the thresholds,
       BWregionprops(i).Area = 0;                                                     % then set the area to zero for regions.                        
    end
end

for i = 1:numel(BWregionprops)
    bbox = BWregionprops(i).BoundingBox;
    topEdge = bbox(2);
    bottomEdge = topEdge + bbox(4) - 1;

    % Remove regions touching the top or bottom
    if topEdge <= 1 || bottomEdge >= imHeight
        BWregionprops(i).Area = 0;
    end
end

BWmask = ismember(labelmatrix(bwconncomp(BWim)), find([BWregionprops.Area] > 0)); % Generate a binary image with filtered objects

BWfill = imfill(BWmask, 'holes');

imThresh          = imX; % Create thresholded image based on points passed by BWmask and original values
imThresh(~BWfill) = 0;
    
end