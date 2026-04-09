function circleLocker = circle_alignment(folder)
% circle_alignment:

%% Section 1: Get list of all TIF files in the folder, including subfolders

% Get list of all TIF files in the folder, including subfolders
imFiles = fullfile(folder, '**', '*.tif');
full    = dir(imFiles);

%% Section 2: Initialize Data Structure to store image results

circleLocker(length(full)) = struct('filename', [], ...
                                    'dna', [], 'refProt', [], 'intProt', [], ...
                                    'intProttranslate', [], 'intProtseg', [], ...
                                    'ogCent', [], 'translateCent', [], ...
                                    'refProtThresh', [], 'dnaThresh2x', []);

imInterval = 1:512;
%% Section 3: Main loop; Process Each Image

% Loop through each image in the folder
for i = 1:length(full)
    %% Section 4: Read & collect max projections from the current image

    % Read the current image
    source = fullfile(full(i).folder, full(i).name);
    data   = three_channel_stacks(source);
    
    % Store the output variables individually for eaiser referencing
    dna     = data.dna;    
    refProt = data.refProt;
    intProt = data.refProt;

     %% Section 5: Segment & Threshold DNA
    
    dnaSeg    = segment_dna_2d(dna); % Perform segmentation on the DNA image using a custom function, 'segment_dna_2d'
    dnaBinSeg = dnaSeg > 0;          % Create a binary mask where dnaThresh is greater than zero (thresholded regions)

    dnaRegionprops = regionprops(dnaBinSeg, 'Area');                                   % Compute region properties (Area) for connected components in the binary mask
    validRegions   = [dnaRegionprops.Area] >= 1000;                                       % Identify valid regions where the area is at least 1000 pixels
    dnaMask        = ismember(labelmatrix(bwconncomp(dnaBinSeg)), find(validRegions)); % Create a mask that includes only valid connected components

    dnaThresh2x = zeros(512); % Initialize another blank 512x512 matrix for thresholded DNA values
    
    % Iterate through image interval
    for x = imInterval
        for y = imInterval
            if dnaMask(x, y) > 0                 % If the dnaMask at (x, y) is part of a valid region,
               dnaThresh2x(x, y) = dnaSeg(x, y); % then retain dnaThresh value
            end
        end
    end
    
    %% Section 6: Compute Image Contrast & Adjust Segmentation Sensitivity
    
    contrast            = std(double(refProt(:))) / mean(double(refProt(:)));      % Compute image contrast metric
    adaptiveSensitivity = max(0.3, 0.5 - 0.2 * contrast);                          % Adjust sensitivity dynamically (tweak scaling factor as needed)
    refProtThresh       = segment_singlecell_circle(refProt, adaptiveSensitivity); % Call the segmentation function with desired sensitivity

    %% Section 7: Create Binary Mask & Extract Centroid
    
    refProtMask = refProtThresh > 0; % Convert to binary mask

    refProtRegionprops = regionprops(refProtMask, 'Centroid', 'Area'); % Find the centroid of the current image

    if isempty(refProtRegionprops)                                         % If no regions were found,
       warning('No refProt regions found in %s. Skipping.', full(i).name); % then display message
       continue;                                                         % and skip this image.
    end
    
    % Find the centroid of the largest region just in case there are multiple.
    [~, largestIdx] = max([refProtRegionprops.Area]); 
    centroids       = refProtRegionprops(largestIdx).Centroid;
        
    %% Section 8: Define Reference Centroid

    % If this is the first image, set the reference centroid
    if i == 1
        refCent = [256, 256];
    end
        
    %% Section 9: Extract IntProt Signal in RefProt Region

    intProtseg = zeros(512, 512); % Masked IntProt signal in refProt region

    % Loop through all pixels to apply masking based on refProt region
    for x = imInterval
        for y = imInterval
            if refProtMask(x,y) == 1            % if refProt mask is true
               intProtseg(x,y)  = intProt(x,y); % then pass intProt for that point
            end
        end
    end
        
    %% Section 10: Align IntProt Image to Reference Centroid
    
    translateVec = refCent - centroids; % Calculate translation vector between the reference centroid and the current image centroid
    
    intProttranslate = imtranslate(intProtseg, translateVec); % Apply the translation to align the IntProt image to the reference centroid

    %% Section 11: Compute & Store Centroids for IntProt Before & After Translation

    intProtmask                = intProtseg > 0;                                % Binary mask for IntProt
    intProttranslateMask       = intProttranslate > 0;                          % Binary mask after translation
    intProtcentStruct          = regionprops(intProtmask, 'Centroid');          % Find original centroid of IntProt
    intProtCent                = intProtcentStruct.Centroid;                    % Store in variable
    intProttranslateCentStruct = regionprops(intProttranslateMask, 'Centroid'); % Find translated centroid
    intProttranslateCent       = intProttranslateCentStruct.Centroid;           % Store in variable

    % Print the old centroid & new centroid
    fprintf('Image: %s\n', full(i).name);
    fprintf('Original IntProt Centroid: (%.2f, %.2f)\n', intProtCent);
    fprintf('Translated IntProt Centroid: (%.2f, %.2f)\n', intProttranslateCent);
        
    %% Section 12: Store Original & Processed Data

    circleLocker(i).filename  = full(i).name;

    circleLocker(i).intProt = intProt;
    circleLocker(i).refProt = refProt;
    circleLocker(i).dna     = dna;

    circleLocker(i).intProttranslate = intProttranslate;
    circleLocker(i).intProtseg       = intProtseg;

    circleLocker(i).ogCent        = intProtCent;
    circleLocker(i).translateCent = intProttranslateCent;

    circleLocker(i).refProtThresh = refProtThresh;
    circleLocker(i).dnaThresh2x   = dnaThresh2x;

end % end main loop

%% Section 13: Save Results

save('centlocker.mat', 'circleLocker'); % Save the locker structure to a MAT file for later use

end % end function

%% Section 14: Image Segmentation Helper Function

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

