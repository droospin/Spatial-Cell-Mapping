function bowLocker = bow_alignment(folder)
% bow_alignment:

%% Section 1: Get list of all TIF files in the folder, including subfolders

% Get list of all TIF files in the folder, including subfolders
imFiles = fullfile(folder, '**', '*.tif');
full    = dir(imFiles);

%% Section 2: Initialize Data Structure to store image results

bowLocker(length(full)) = struct('filename', [], ...
                                 'dna', [], 'refProt', [], 'intProt', [], ...
                                 'dnaThresh', [], 'refProtThresh', [], ...
                                 'intProtseg', [], 'normIntProt', [], ...
                                 'refCent', []);

imInterval = 1:512; 
%% Section 3: Read & collect max projection from DNA channel from all images in folder

dnaStruct(length(full)) = struct('dna', zeros(512));

for i = 1:length(full)
    source           = fullfile(full(i).folder, full(i).name);
    data             = three_channel_stacks(source);
    dnaStruct(i).dna = data.dna;
    disp(['Processing sample ', num2str(i), ': ', full(i).name]);
end

%% Section 4: Main loop; Process Each Image

for i = 1:length(full)
    %% Section 5: Read & collect max projections from the current image
    
    dnaCurr = dnaStruct(i).dna; % Retrieve pre-loaded DNA data for the specific iteration

    % Skip iteration if max intensity is below 4000
    if max(dnaCurr(:)) < 4000
        disp(['Skipping sample ', num2str(i), ' due to low DNA intensity.']);
        continue;  % Skip to the next iteration
    end
    
    % Loop through each image in the folder
    source = fullfile(full(i).folder, full(i).name);
    data   = three_channel_stacks(source);

    refProt  = data.refProt; 
    intProt  = data.refProt;

    %% Section 6: Segment & Threshold DNA

    dnaThresh = zeros(512); % Initialize a blank 512x512 matrix to store high-intensity DNA values

    % Iterate through the given image interval
    for x = imInterval
        for y = imInterval
            if dnaCurr(x, y) > 3000             % If the DNA intensity at (x, y) is greater than 3000,
               dnaThresh(x, y) = dnaCurr(x, y); % then retain it in realDna
            end
        end
    end

    dnaThresh = segment_dna_2d(dnaThresh);

    binaryDNA = dnaThresh > 0;
    dnaRegionprops = regionprops(binaryDNA, 'Area');
    for n = 1:numel(dnaRegionprops)
        if dnaRegionprops(n).Area < 1000
           dnaRegionprops(n).Area = 0;
        end
    end

    dnaThreshMask = ismember(labelmatrix(bwconncomp(binaryDNA)), find([dnaRegionprops.Area] > 0));
    dnaThresh(~dnaThreshMask) = 0;

    %% Section 7: Run rotate_and_crop helper function

    [dnaThresh, refProt, intProt] = rotate_and_crop(0, 50, dnaThresh, refProt, intProt); 

    %% Section 10: Compute Image Contrast & Adjust Segmentation Sensitivity
    
    contrast            = std(double(refProt(:))) / mean(double(refProt(:)));   % Compute image contrast metric
    adaptiveSensitivity = max(0.1, 0.3 - 0.2 * contrast);                       % Adjust sensitivity dynamically (tweak scaling factor as needed)
    refProtThresh       = segment_singlecell_bow(refProt, adaptiveSensitivity); % Call the segmentation function with desired sensitivity

    %% Section 11: Find RefProt Centroid

    binaryRefProt      = refProtThresh > 0; % Create a binary mask where refProtThresh is greater than 0
    
    refProtRegionprops = regionprops(binaryRefProt, 'Centroid');      % Compute the centroid of connected components in the binary mask
    disp(refProtRegionprops)                                          % Display the properties of the detected regions
    disp(['Processing sample ', num2str(i), ': ', full(i).name]);     % Print the current iteration number
    refProtRegionprops.Centroid = round(refProtRegionprops.Centroid); % Round the centroid coordinates to the nearest integer
    
    %%

    imshowpair(dnaThresh > 0, refProtThresh > 0) % Display the binary thresholded DNA image
    hold on;
    plot(refProtRegionprops.Centroid(1), refProtRegionprops.Centroid(2), 'rx', 'MarkerSize', 10) % Plot the DNA centroid as a red 'X' marker
    hold off;
    pause(0.5); % Pause for 0.5 seconds to allow visualization
    
    %% Section 11: Extract IntProt Signal in RefProt Region

    refProtMask = refProtThresh > 0;
    intProtseg  = zeros(512, 512); % Masked IntProt signal in refProt region

    % Loop through all pixels to apply masking based on refProt region
    for x = imInterval
        for y = imInterval
            if refProtMask(x,y) == 1
               intProtseg(x,y)   = intProt(x,y) - 4000;
            end
        end
    end

    %% Section 12: Store all necessary & important values in a structure

    intProtmax  = max(intProtseg(:));      % find intProt max intensity to normalize
    intProtNorm = intProtseg / intProtmax; % normalize intProt based off max intensity
    
    bowLocker(i).filename = full(i).name;

    bowLocker(i).intProt = intProt;
    bowLocker(i).refProt = refProt;
    bowLocker(i).dna     = dnaThresh;

    bowLocker(i).intProtseg    = intProtseg;
    bowLocker(i).refProtThresh = refProtThresh;

    bowLocker(i).intProtNorm = intProtNorm;
    bowLocker(i).refCent     = refProtRegionprops.Centroid; 

    clearvars('-except', "bowLocker", "full", "imageFiles", "folder", "dnaStruct", "imInterval");
    close all;

end % main loop
end % full function

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

