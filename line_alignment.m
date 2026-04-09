function lineLocker = line_alignment(folder)
% line_alignment:

%% Section 1: Get list of all TIF files in the folder, including subfolders

% Get list of all TIF files in the folder, including subfolders
imFiles = fullfile(folder, '**', '*.tif');
full    = dir(imFiles);

%% Section 2: Initialize Data Structure to store image results

lineLocker(length(full)) = struct('filename', [], ...
                                  'dna', [], 'refProt', [], 'intProt', [], ...
                                  'dnaThresh2x', [], 'refProtThresh', [], ...
                                  'intProtseg', [], 'normIntProt', [], ...
                                  'dnaCent', []);

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
    
    if ismember(i, [1 2 3 4 5 6 7 8 9])
        % Store the output variables individually for easier referencing
        refProt  = data.refProt;
        intProt  = data.refProt;
    else
        refProt = data.refProt;
        intProt = data.refProt;
    end

    %% Section 6: Segment & Threshold DNA

    dnaReal = zeros(512); % Initialize a blank 512x512 matrix to store high-intensity DNA values

    % Iterate through the given image interval
    for x = imInterval
        for y = imInterval
            if dnaCurr(x, y) > 3000           % If the DNA intensity at (x, y) is greater than 3000,
               dnaReal(x, y) = dnaCurr(x, y); % then retain it in realDna
            end
        end
    end

    dnaThresh    = segment_dna_2d(dnaReal); % Perform segmentation on the DNA image using a custom function, 'segment_dna_2d'
    dnaBinThresh = dnaThresh > 0;           % Create a binary mask where dnaThresh is greater than zero (thresholded regions)
    
    dnaRegionprops = regionprops(dnaBinThresh, 'Area');                                   % Compute region properties (Area) for connected components in the binary mask
    validRegions   = [dnaRegionprops.Area] >= 1000;                                       % Identify valid regions where the area is at least 1000 pixels
    dnaMask        = ismember(labelmatrix(bwconncomp(dnaBinThresh)), find(validRegions)); % Create a mask that includes only valid connected components

    dnaThresh2x = zeros(512); % Initialize another blank 512x512 matrix for thresholded DNA values
    
    % Iterate through image interval
    for x = imInterval
        for y = imInterval
            if dnaMask(x, y) > 0                    % If the dnaMask at (x, y) is part of a valid region,
               dnaThresh2x(x, y) = dnaThresh(x, y); % then retain dnaThresh value
            end
        end
    end
    
    %% Section 7: Run rotate_and_crop helper function

    [dnaReal, dnaThresh2x, refProt, intProt] = rotate_and_crop(0, 50, dnaReal, dnaThresh2x, refProt, intProt);  % -5 for 15 & 20, 
    
    %% Section 8: Find dna Centroid

    binaryDNA = dnaThresh2x > 0;                                  % Create a binary mask where dnaThresh2x is greater than 0
    dnaRegionprops = regionprops(binaryDNA, 'Centroid');          % Compute the centroid of connected components in the binary mask
    disp(dnaRegionprops)                                          % Display the properties of the detected regions
    disp(['Processing sample ', num2str(i), ': ', full(i).name]); % Print the current iteration number
    dnaRegionprops.Centroid = round(dnaRegionprops.Centroid);     % Round the centroid coordinates to the nearest integer

    %% Section 9: Flipping cells so that retrrefProtg edge of cell is on left side of the image
    
    if dnaRegionprops.Centroid(1) > 256 || (i == 6)
                                         % Check if the x-coordinate of the DNA centroid is greater than 256            
                                         % flip those meeting the qualifiications

                                         % FOR 5 µm:  && ~(i == 1 || i == 12) || (i == 4 || i == 9 || i == 10)
                                         % FOR 10 µm: || (i == 6)
                                         % FOR 15 µm: && ~(i == 4 || i == 11)


       disp(['Flipping at iteration ', num2str(i)]); % Display a message indicating that flipping is occurring
       dnaReal     = flip(dnaReal, 2);               % Flip the images horizontally (along the second dimension)
       dnaThresh2x = flip(dnaThresh2x, 2);
       refProt     = flip(refProt, 2);
       intProt     = flip(intProt, 2);
    
       dnaRegionprops.Centroid(1) = 512 - dnaRegionprops.Centroid(1); % Adjust the centroid's x-coordinate after flipping
    end
    
    imshow(dnaThresh2x > 0) % Display the binary thresholded DNA image
    hold on;
    plot(dnaRegionprops.Centroid(1), dnaRegionprops.Centroid(2), 'rx', 'MarkerSize', 10) % Plot the DNA centroid as a red 'X' marker
    hold off;
    pause(0.5); % Pause for 0.5 seconds to allow visualization

    %% Section 10: Compute Image Contrast & Adjust Segmentation Sensitivity
    
    contrast            = std(double(refProt(:))) / mean(double(refProt(:)));    % Compute image contrast metric
    adaptiveSensitivity = max(0.1, 0.3 - 0.2 * contrast);                        % Adjust sensitivity dynamically (tweak scaling factor as needed)
    refProtThresh       = segment_singlecell_line(refProt, adaptiveSensitivity); % Call the segmentation function with desired sensitivity

    %% Section 11: Extract IntProt Signal in RefProt Region

    refProtMask = refProtThresh > 0;
    intProtseg   = zeros(512, 512); % Masked IntProt signal in refProt region

    % Loop through all pixels to apply masking based on refProt region
    for x = imInterval
        for y = imInterval
            if refProtMask(x,y) == 1
               intProtseg(x,y)   = intProt(x,y) - 4000;
            end
        end
    end

    %% Section 12: Store all necessary & important values in a structure

    intProtmax  = max(intProtseg(:));   % find intProt max intensity to normalize
    intProtNorm = intProtseg / intProtmax; % normalize intProt based off max intensity
    
    lineLocker(i).filename = full(i).name;

    lineLocker(i).intProt  = intProt;
    lineLocker(i).refProt = refProt;
    lineLocker(i).dna   = dnaReal;

    lineLocker(i).intProtseg     = intProtseg;
    lineLocker(i).refProtThresh = refProtThresh;
    lineLocker(i).dnaThresh2x = dnaThresh2x;

    lineLocker(i).intProtNorm = intProtNorm;
    lineLocker(i).dnaCent  = dnaRegionprops.Centroid; 

    clearvars('-except', "lineLocker", "full", "imageFiles", "folder", "dnaStruct", "imInterval");
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

