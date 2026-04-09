function [lineLocker] = line_3d(folder)

%% Section 1: Get list of all TIF files in the folder, including subfolders

% Define specific subfolder names
subfolderNames = {'5um', '10um', '15um', '20um'};

% Initialize
full         = []; % will be a struct array of file info
folderLabels = {}; % to track folder for each file

% Loop through each target subfolder
for k = 1:length(subfolderNames)
    subfolderPath = fullfile(folder, subfolderNames{k});
    imFiles       = dir(fullfile(subfolderPath, '**', '*.tif'));
    % Assign full path and folder label
    for j = 1:length(imFiles)
        folderLabels{end+1,1} = subfolderNames{k};
    end
    
    full = [full; imFiles];
end

numCells     = numel(full);
folderLabels = string(folderLabels);  % convert to string for indexing

% --- [NEW] Track relative index within each subfolder ---
folderCounts  = containers.Map({'5um', '10um', '15um', '20um'}, [0, 0, 0, 0]);
relativeIndex = zeros(numCells, 1);

for i = 1:numCells
    label               = folderLabels(i);
    folderCounts(label) = folderCounts(label) + 1;
    relativeIndex(i)    = folderCounts(label);
end

%% Section 2: Establish constant variables for main loop and function in general

% Create custom 'hot' color channel
numColors    = 256;                                                              % Define number of color levels
origColormap = hot(numColors);                                                   % Get the default 'hot' colormap
cutoff       = 0.9;                                                              % Define a cutoff threshold (e.g., 80% of the colormap)
newMaxIdx    = round(cutoff * numColors); 
customHot    = origColormap(1:newMaxIdx, :);                                     % Stretch and cut off colormap
customHot    = [customHot; repmat(customHot(end, :), numColors - newMaxIdx, 1)]; % Extend the last color to fill up to 1


imInterval     = 1:512;                     % image size for x & y dimensions for looping 
numBins        = 100;                       % or whatever number of bins you want
fixedDistances = linspace(-1, 1, numBins);  % creating spacing vector based on desired 'numBins' from -1 to 1
Ntarget        = 50;                        % or your desired z-slices

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
cellDepths = nan(numCells, 1);
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
allRefMean = nan(numCells, Ntarget, numBins); % one layer per cell
allIntMean = nan(numCells, Ntarget, numBins);

%% -----------------------------
%% 2) Second pass: process each cell and interpolate onto the common z‐grid
%% -----------------------------
for i = 1:numCells
    %% ---- read & segment exactly as before ----
    source = fullfile(full(i).folder, full(i).name);

    disp(full(i).name) % display file name
    
    data               = three_channel_stacks(source);
    [xDim, yDim, zDim] = size(data.stackDNA);

    label = string(folderLabels(i));

    switch label
        case {"5um","10um"}
            dna        = data.stackDNA;
            refProt    = data.stackRefProt;
            intProt    = data.stackIntProt;
            maxDNA     = data.dna;
            maxRefProt = data.refProt;
            rotateAngle = 0;
    
        case {"15um","20um"}
            dna        = data.stackIntProt;
            refProt    = data.stackDNA;
            intProt    = data.stackRefProt;
            maxDNA     = data.intProt;
            maxRefProt = data.dna;
            rotateAngle = -5;
    
        otherwise
            error("Unknown folder logic for %s", label);
    end

    %% Section 7: Run rotate_and_crop helper function

    [dna, refProt, intProt, maxDNA, maxRefProt] = rotate_and_crop(rotateAngle, 50, dna, refProt, intProt, maxDNA, maxRefProt);  % -5 for 15 & 20, 
    %% Section 6: Segment & Threshold DNA

    dnaReal = zeros(512); % Initialize a blank 512x512 matrix to store high-intensity DNA values

    % Iterate through the given image interval
    for x = imInterval
        for y = imInterval
            if maxDNA(x, y) > 3000           % If the DNA intensity at (x, y) is greater than 3000,
               dnaReal(x, y) = maxDNA(x, y); % then retain it in realDna
            end
        end
    end

    maxDNAthresh = segment_dna_2d(dnaReal); % Perform segmentation on the DNA image using a custom function, 'segment_dna_2d'
    dnaBinThresh = maxDNAthresh > 0;           % Create a binary mask where dnaThresh is greater than zero (thresholded regions)
    
    dnaRegionprops = regionprops(dnaBinThresh, 'Area');                                   % Compute region properties (Area) for connected components in the binary mask
    validRegions   = [dnaRegionprops.Area] >= 1000;                                       % Identify valid regions where the area is at least 1000 pixels
    dnaMask        = ismember(labelmatrix(bwconncomp(dnaBinThresh)), find(validRegions)); % Create a mask that includes only valid connected components

    dnaThresh2x = zeros(512); % Initialize another blank 512x512 matrix for thresholded DNA values
    
    % Iterate through image interval
    for x = imInterval
        for y = imInterval
            if dnaMask(x, y) > 0                       % If the dnaMask at (x, y) is part of a valid region,
               dnaThresh2x(x, y) = maxDNAthresh(x, y); % then retain dnaThresh value
            end
        end
    end

    imshow(dnaThresh2x, []); pause(0.5);
    %% Section 8: Find dna Centroid

    binaryDNA      = dnaThresh2x > 0;                             % Create a binary mask where dnaThresh2x is greater than 0
    dnaRegionprops = regionprops(binaryDNA, 'Centroid', 'Area');          % Compute the centroid of connected components in the binary mask
    disp(dnaRegionprops)                                          % Display the properties of the detected regions
    disp(['Processing sample ', num2str(i), ': ', full(i).name]); % Print the current iteration number

    [~, largestIdx]         = max([dnaRegionprops.Area]); 
    dnaRegionprops.Centroid = dnaRegionprops(largestIdx).Centroid;
    dnaRegionprops.Centroid = round(dnaRegionprops.Centroid);     % Round the centroid coordinates to the nearest integer

     %% Section 9: Flipping cells so that retrrefProtg edge of cell is on left side of the image
    
    flipThis = false;
    % Get x-coordinate of centroid
    xCentroid = dnaRegionprops.Centroid(1);

    relIdx = relativeIndex(i);  % <-- [NEW] subfolder-relative index
    
    % Folder-specific flipping logic
    switch folderLabels(i)
        case '5um'
            if (xCentroid > 256 && ~(relIdx == 1 || relIdx == 12)) || ismember(relIdx, [4, 9, 10])
               flipThis = true;
            end
        case '10um'
            if xCentroid > 256 || relIdx == 6
               flipThis = true;
            end
        case '15um'
            if xCentroid > 256 && ~(relIdx == 4 || relIdx == 11)
               flipThis = true;
            end
        case '20um'
            if xCentroid > 256
               flipThis = true;
            end
    end
    
    % Perform flipping if condition is met
    if flipThis
        disp(['Flipping at iteration ', num2str(i)]); % Display a message
    
        dnaReal     = flip(dnaReal, 2);
        dnaThresh2x = flip(dnaThresh2x, 2);
        dna         = flip(dna, 2);
        maxRefProt  = flip(maxRefProt, 2);
        refProt     = flip(refProt, 2);
        intProt     = flip(intProt, 2);
    
        % Adjust centroid after flipping
        dnaRegionprops.Centroid(1) = 512 - xCentroid;
    end
    
    imshowpair(dnaThresh2x, intProt(:, :, 11)); pause(0.5); 

    dnaCent = dnaRegionprops.Centroid;

    %%

    dnaSeg = zeros(size(dna));
    for x = 1:xDim
        for y = 1:yDim
            for z = 1:zDim
                if dnaThresh2x(x, y) > 0
                   dnaSeg(x, y, z) = dnaThresh2x(x, y);
                end
            end
        end
    end

    for x = 1:xDim
        for y = 1:yDim
            for z = 1:zDim
                if dnaSeg(x, y, z) < 3500
                   dnaSeg(x, y, z) = 0;
                end
            end
        end
    end

    dnaBinSeg      = dnaSeg > 0;
    dnaRegionprops = regionprops(dnaBinSeg, 'Area');
    validRegions   = [dnaRegionprops.Area] >= 2000;
    dnaMask        = ismember(labelmatrix(bwconncomp(dnaBinSeg)), find(validRegions));

    dnaThresh = zeros(size(dna));

    for x = 1:xDim
        for y = 1:yDim
            for z = 1:zDim
                if dnaMask(x, y, z) > 0
                   dnaThresh(x, y, z) = dnaSeg(x, y, z);
                end
            end
        end
    end
    
    %% Section 10: Compute Image Contrast & Adjust Segmentation Sensitivity
    
    contrast            = std(double(maxRefProt(:))) / mean(double(maxRefProt(:))); % Compute image contrast metric
    adaptiveSensitivity = max(0.1, 0.3 - 0.2 * contrast);                           % Adjust sensitivity dynamically (tweak scaling factor as needed)
    maxRefProtThresh    = segment_singlecell_line(maxRefProt, adaptiveSensitivity); % Call the segmentation function with desired sensitivity  % Call the segmentation function with desired sensitivity
    
    imshowpair(dnaThresh2x, maxRefProtThresh) % Display the binary thresholded DNA image
    hold on;
    plot(dnaCent(1), dnaCent(2), 'rx', 'MarkerSize', 10) % Plot the DNA centroid as a red 'X' marker
    hold off;
    pause(1); % Pause for 0.5 seconds to allow visualization

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
                if refProtSeg(x, y, z) < 3000
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
                if intProtSeg(x, y, z) < 3000
                   intProtSeg(x, y, z) = 0;
                end
            end
        end
    end

    imshow(intProtSeg(:, :, 7), []); pause (0.5);

    %%

    refProtMask   = refProtSeg > 0;
    refProtThresh = zeros(size(refProtMask));

    for z = 1:size(refProtSeg, 3)
        sliceMask            = refProtMask(:,:,z);
        cc                   = bwconncomp(sliceMask);
        stats                = regionprops(cc, 'Area');
        validIdx             = find([stats.Area] >= 100);
        cleanMask            = ismember(labelmatrix(cc), validIdx);
        refProtMask(:, :, z) = cleanMask;
    
        slice                  = double(refProtSeg(:,:,z));  % <-- convert to double
        mask                   = double(cleanMask);          % <-- also convert mask
        refProtThresh(:, :, z) = slice .* mask;              % <-- safe to multiply
    end

   %% Compute Euc Dist from nuclear center for all pixels in each slice —

   REFeucDistNucCent = nan(xDim, yDim, zDim);

   for y = imInterval
        for x = imInterval
            for z = 1:size(refProtThresh, 3)
                if refProtThresh(y, x, z) > 0
                   REFeucDistNucCent(y, x, z) = sqrt((x - dnaCent(1))^2 + (y - dnaCent(2))^2);

                   if x < dnaCent(1)
                      REFeucDistNucCent(y, x, z) = - REFeucDistNucCent(y, x, z);
                   end

                end
            end
        end
   end

   imshow(refProtThresh(:, :, 7), []); pause(0.5);
   imshow(REFeucDistNucCent(:, :, 7)); pause(0.5);

    %% -----------------------------
    %% 2a) Normalize radial distances and intensities (unchanged)
    %% -----------------------------

    maxRefEucDist  = max(REFeucDistNucCent(~isnan(REFeucDistNucCent)));

    % initialize
    normRefEucDist = nan(size(REFeucDistNucCent));

    normRefInt = zeros(size(refProtThresh));
    normIntInt = zeros(size(intProtSeg));

    % annotate valid pixels
    validMask = refProtThresh > 0  & ~isnan(REFeucDistNucCent);

    normRefEucDist(validMask) = REFeucDistNucCent(validMask) / maxRefEucDist;

    normRefInt(validMask) = refProtThresh(validMask) / max(refProtThresh(:));
    normIntInt(validMask) = intProtSeg(validMask)    / max(intProtSeg(:));

    sum(normRefEucDist(:) == 0)

    %{
    minRefDist = abs(min(normRefEucDist(:)));
    for y = imInterval
        for x = imInterval
            for z = 1:size(refProt, 3)
                if normRefEucDist(y, x, z) < 0
                   normRefEucDist(y, x, z) = normRefEucDist(y, x, z) / minRefDist;
                end
                if normIntEucDist(y, x, z) < 0
                   normIntEucDist(y, x, z) = normIntEucDist(y, x, z) / minRefDist;
                end
            end
        end
    end
    %}

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
    %% 2d) INTERPOLATE each cell's (Zi × numBins) onto (Ntarget × numBins)
    %% (Option 1: normalize z‐axis).
    %% -----------------------------
    % Build normalized z‐coordinate for this cell:
    % If Zi == 1, we'll just replicate that row across all targets.
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
            % Use linear interpolation; "extrap" set to NaN by default if out of range
            refInterp(:, b) = interp1(t_orig, v_ref, t_target, 'linear', NaN);
            intInterp(:, b) = interp1(t_orig, v_int, t_target, 'linear', NaN);
        end
    end

    % Store into the big 3D accumulator (cell i's layer)
    allRefMean(i,:,:) = refInterp;  
    allIntMean(i,:,:) = intInterp;

        %% Section 7: Store Original & Processed Data

    lineLocker(i).filename = full(i).name;

    lineLocker(i).refProt = refProt;
    lineLocker(i).intProt = intProt;
    lineLocker(i).dna     = dna;

    lineLocker(i).intProtSeg = intProtSeg;
    lineLocker(i).refProtSeg = refProtSeg;

    lineLocker(i).refProttranslate = refProtSeg;
    lineLocker(i).intProttranslate = intProtSeg;

    lineLocker(i).refProtThresh = refProtThresh;
    lineLocker(i).dnaThresh     = dnaThresh;

    lineLocker(i).dnaCent = dnaCent;

    lineLocker(i).REFeucDistCent = REFeucDistNucCent;

    lineLocker(i).normRefEucDist = normRefEucDist;

    lineLocker(i).rawRefSideProfile = refMeanInt_i; % size Zi×numBins
    lineLocker(i).rawIntSideProfile = intMeanInt_i; % size Zi×numBins

    lineLocker(i).allRefMean = allRefMean;
    lineLocker(i).allIntMean = allIntMean;

    lineLocker(i).depth = cellDepths(i);

    lineLocker(i).fixedDistances = fixedDistances;
    lineLocker(i).Ntarget        = Ntarget;
    lineLocker(i).t_target       = linspace(0,1,Ntarget);

    lineLocker(i).thresholds.dna        = 5500;
    lineLocker(i).thresholds.refProt    = 3500;
    lineLocker(i).thresholds.minAreaDNA = 2000;
    lineLocker(i).thresholds.minAreaRef = 100;

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

    lineLocker(1).avgRefMean = avgRefMean;
    lineLocker(1).avgIntMean = avgIntMean;
    
    %% -----------------------------
    %% 4) Plot the final “mean across cells” images
    %% -----------------------------
    
    figure;
    imagesc(fixedDistances, linspace(1, Ntarget, Ntarget), avgRefMean);
    set(gca, 'YDir','normal');
    %xlabel('Normalized Distance (bin)');
    %ylabel('Normalized Slice Index');
    %title('⟨Reference Channel⟩ – Averaged Over All Cells');
    colormap(customHot);
    %colorbar;
    axis off;
    
    figure;
    imagesc(fixedDistances, linspace(1, Ntarget, Ntarget), avgIntMean);
    set(gca, 'YDir','normal');
    %xlabel('Normalized Distance (bin)');
    %ylabel('Normalized Slice Index');
    %title('⟨Interest Channel⟩ – Averaged Over All Cells');
    colormap(customHot);
    %colorbar;
    axis off;

    % end main loop
end

%% Section 13: Helper Function to rotate images based on orientation to x-axis in raw images

function varargout = rotate_and_crop(angle, padSize, varargin)

numImages = numel(varargin);
varargout = cell(1, numImages);

for n = 1:numImages

    im = varargin{n};

    if ndims(im) == 2
        % ---- 2D case ----
        padIm    = padarray(im, [padSize padSize], 'replicate', 'both');
        rotIm    = imrotate(padIm, angle, 'bilinear', 'crop');
        outIm    = rotIm(padSize+1:end-padSize, padSize+1:end-padSize);

    elseif ndims(im) == 3
        % ---- 3D case: rotate each Z slice ----
        [H, W, Z] = size(im);
        outIm = zeros(H, W, Z, 'like', im);

        for z = 1:Z
            padSlice = padarray(im(:,:,z), [padSize padSize], 'replicate', 'both');
            rotSlice = imrotate(padSlice, angle, 'bilinear', 'crop');
            outIm(:,:,z) = rotSlice(padSize+1:end-padSize, padSize+1:end-padSize);
        end

    else
        error('rotate_and_crop only supports 2D or 3D arrays.');
    end

    varargout{n} = outIm;

end
end


%% Section 14: Helper function to threshold cell mask channel

function imThresh = segment_singlecell_line(imX, sensitivity)
% Segment image using auto-generated code from Image Segmenter app
% [BW,MASKEDIMAGE] = segmentImage(X) segments image X using auto-generated
% code from the Image Segmenter app. The final segmentation is returned in
% BW, and a masked image is returned in MASKEDIMAGE.
    
% Partly auto-generated by imageSegmenter app on 17-Feb-2025
%----------------------------------------------------

imX = imadjust(imX); % Adjust data to span data range
    
BWim = imbinarize(im2gray(imX), 'adaptive', 'Sensitivity', sensitivity, 'ForegroundPolarity', 'dark'); % Threshold image with adaptive threshold
    
BWfill = imfill(BWim, 'holes'); % Fill holes

BWregionprops = regionprops(BWfill, 'Area', 'BoundingBox'); % Compute properties of connected components in the binary image BWfill

minSizeThresh = 5000;  % Define minimum
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

BWmask = ismember(labelmatrix(bwconncomp(BWfill)), find([BWregionprops.Area] > 0)); % Generate a binary image with filtered objects

imThresh          = imX; % Create thresholded image based on points passed by BWmask and original values
imThresh(~BWmask) = 0;
    
end