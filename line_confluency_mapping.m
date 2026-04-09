function lineOut = line_confluency_mapping(folder)
% This function analyzes multi-channel fluorescence microscopy images to map
% cell confluency and protein localization. It processes a stack of images by:
% - Thresholding and masking reference and intensity markers
% - Segmenting nuclei and identifying cell boundaries using watershed
% - Creating pattern edge masks
% - Computing distance maps within each cell
% - Classifying protein intensities as junctional, edge-associated, or internal
% The results are returned as a structured array for further analysis.

%% Section 1: Get list of all TIF files in the folder, including subfolders

% Get list of all TIF files in the folder, including subfolders
imFiles = fullfile(folder, '**', '*.tif');
full    = dir(imFiles);

%% Main Loop

for i = 1:length(full)

    disp(i);
    disp(full(i).name);
    
%% Image reading
    source = fullfile(full(i).folder, full(i).name);
    data   = four_channel_stacks(source);
    
    refProt  = data.refProt;
    intProt1 = data.intProt1;
    intProt2 = data.intProt2;
    dna      = data.dna;

%% general thresholding

   % refProt(refProt < 4500) = 0;
   % intProt1(intProt1 < 5500) = 0;

%% dna thresholding

    dnaThresh     = segment_dna_2d(dna);
    dnaMask1      = dnaThresh > 0;
    minThreshSize = 500;

    dnaRegionprops = regionprops(dnaMask1, 'Area');

    for k = 1:numel(dnaRegionprops)
            if dnaRegionprops(k).Area < minThreshSize % If an area is less than the 'minThreshSize,'
               dnaRegionprops(k).Area = 0;            % then delete areas smaller than 'minSizeThresh'
            end
    end
           
    % remove areas smaller than 'minSizeThresh' in our image
    dnaMask2 = ismember(labelmatrix(bwconncomp(dnaMask1)), find([dnaRegionprops.Area] > 0));
    
    dnaThresh2x = dnaThresh;       % copy full image
    dnaThresh2x(dnaMask2 == 0) = 0;  % zero out values outside the mask
    
    dnaMask3 = dnaThresh2x > 0;

%%

    imshowpair(intProt1, refProt);
    title('Raw intProt1 (green) & refProt (magenta)');
    pause(1);

%% Create Pattern Mask

    pattMask1     = refProt > 3500;
    minThreshSize = 10000;
    
    pattRegionprops = regionprops(pattMask1, 'Area');
    
    for k = 1:numel(pattRegionprops)
            if pattRegionprops(k).Area < minThreshSize % If an area is less than the 'minThreshSize,'
               pattRegionprops(k).Area = 0;            % then delete areas smaller than 'minSizeThresh'
            end
    end
           
    % remove areas smaller than 'minSizeThresh' in our image
    pattMask2 = ismember(labelmatrix(bwconncomp(pattMask1)), find([pattRegionprops.Area] > 0));
    
    pattMask = imfill(pattMask2, "holes");
    
    imshow(pattMask)
    title('Pattern Mask (logical)');
    pause(1);

%% Find pattern edges (top and bottom)

    % Ensure pattMask is logical
    pattMask = logical(pattMask);
    
    % Initialize blank edge masks
    topEdge1    = false(size(pattMask));
    bottomEdge1 = false(size(pattMask));
    
    % Loop over columns
    for col = 1:size(pattMask, 2)
        colData    = pattMask(:, col);
        rowIndices = find(colData);
        if ~isempty(rowIndices)
            topEdge1(rowIndices(1), col)      = true;           % Topmost pixel in column
            bottomEdge1(rowIndices(end), col) = true;           % Bottommost pixel in column
        end
    end
    
    % Dilate line to ensure full edge of pattern isn't counted later on
    se         = strel('disk', 14, 0);
    topEdge    = imdilate(topEdge1, se);
    bottomEdge = imdilate(bottomEdge1, se);
    
    pattEdge = double(topEdge | bottomEdge);  % cleaner
    
    imshowpair(pattEdge, refProt);
    title('Pattern Edge (green) & refProt (magenta)');
    pause(1);

%% refProt thresholding 

    % Adjust data to span data range.
    refProtAdj = imadjust(refProt);
    
    % Threshold image with global threshold
    BWrefProt = imbinarize(im2gray(refProtAdj));
    
    % Create masked image.
    refProtThresh = refProtAdj;
    refProtThresh(~BWrefProt) = 0;
    
    imshow(refProtThresh)
    
    refProtMask = refProtThresh > 0;
    
    minThreshSize = 500;  % 10000 for line
    
    refProtRegionprops = regionprops(refProtMask, 'Area');
    
    for k = 1:numel(refProtRegionprops)
            if refProtRegionprops(k).Area < minThreshSize % If an area is less than the 'minThreshSize,'
               refProtRegionprops(k).Area = 0;            % then delete areas smaller than 'minSizeThresh'
            end
    end
               
    % remove areas smaller than 'minSizeThresh' in our image
    refProtMask2 = ismember(labelmatrix(bwconncomp(refProtMask)), find([refProtRegionprops.Area] > 0));
    
    refProtThresh(~refProtMask2) = 0;
    imshow(refProtThresh, []);
    title('refProtThresh');
    pause(1);

    close all;

%% Use remove_polygon helper function on refProtThresh

    refProtThreshEdge = remove_polygon(refProtThresh);

%% Watershed cell segmentation

    refWater = watershed(imgaussfilt(refProtThreshEdge, 3)); % 3 is good for lines
    refMask = mergeRegionsGUI(refWater, refProt, dnaMask3);


%% Dilate mask to merge regions
    
    se = strel('disk', 1);  % Adjust size as needed
    cellMask = zeros(size(refMask), 'like', refMask);
    
    labels = unique(refMask);
    labels(labels == 0) = [];  % Remove background
    
    for j = 1:length(labels)
        label = labels(j);
        singleRegion = refMask == label;
        
        % Dilate this region
        dilatedRegion = imdilate(singleRegion, se);
        
        % Remove overlaps with already-filled pixels
        dilatedRegion(cellMask > 0) = 0;
        
        % Assign the label
        cellMask(dilatedRegion) = label;
    end
    
    % Optional: Visualize the dilated labeled mask
    imshow(label2rgb(cellMask, 'jet', 'k'))
    title('Individually Dilated CellMask');
    pause(1);

%% refProt & intProt1 segmenting for pattern

    refPatt = refProt;
    refPatt(pattMask == 0) = 0;
    intPatt = intProt1;
    intPatt(pattMask == 0) = 0;
    
    imshowpair(intPatt, refPatt)
    title('intProt1 Patterned (green) & refProt Patterned (magenta)');
    pause(1);

%% intProt1 segmenting dist vs. intensity per cell

    labels = unique(cellMask);
    labels(labels == 0) = [];  % Remove background

    % Remove small regions less than 4000 pixel area
    for n = 1:numel(labels)
        thisLabel = labels(n);
        area      = sum(cellMask(:) == thisLabel);
        if area < 2000
            cellMask(cellMask == thisLabel) = 0;  % Remove small region
        end
    end
    
    cellPixelIdx = cell(size(labels));
    for n = 1:numel(labels)
        cellPixelIdx{n} = find(cellMask == labels(n));  % indices for each cell
    end

    dists = cell(size(cellPixelIdx));
    for n = 1:numel(cellPixelIdx)
        thisMask = false(size(cellMask));
        thisMask(cellPixelIdx{n}) = true;
    
        edge          = bwperim(thisMask);
        distMap       = bwdist(edge);             % Distance from edge inward
        dists{n}      = distMap(thisMask);        % Extract distances only within the mask
        roundDists{n} = round(dists{n});          % Round those distances
        normDists{n}  = dists{n} / max(dists{n}); % Normalize
    end

    % Initialize output images
    normDistIm  = zeros(size(cellMask));
    roundDistIm = zeros(size(cellMask));
    
    % Populate both images from per-cell distance values
    for n = 1:numel(cellPixelIdx)
        normDistIm(cellPixelIdx{n})  = normDists{n};
        roundDistIm(cellPixelIdx{n}) = roundDists{n};
    end
    
    % Show both images side by side
    figure;
    subplot(1, 2, 1);
    imshow(normDistIm, []);
    title('intProt1 Normalized Distance (0–1)');
    
    subplot(1, 2, 2);
    imshow(roundDistIm, []);
    title('intProt1 Rounded Distance');
    
    pause(1);

%% Create bins for intProt1 (only traced cells)

    intCells = double(intPatt) / max(double(intPatt(:)));
    intCells(normDistIm == 0) = 0; 

    figure; imshow(intCells, []);
    title('intProt1 in isolated cells');
    pause(2);

%  for edge of patt (w/in 10 pixels of pattEdge?)

    normIntProt1pattEdge = intCells;
    normIntProt1pattEdge(pattEdge == 0 | refProtThreshEdge == 0) = 0;

    figure; imshow(normIntProt1pattEdge, []);
    title('intProt1 at the Pattern Edge');
    pause(2);
    
%  junctional places (edge of cell & not pattEdge)

    normIntProt1junc = intCells;
    normIntProt1junc(pattEdge == 1 | roundDistIm >= 10 | refProtThreshEdge == 0) = 0;

    figure; imshow(normIntProt1junc, [])
    title('Junctional intProt1 (not at the edge of the pattern)');
    pause(2);

%  for edge of cell (both edge of patt & junctional)

    normIntProt1edge = intCells;
    normIntProt1edge(roundDistIm >= 10 | refProtThreshEdge == 0) = 0;

    figure; imshow(normIntProt1edge, []);
    title('intProt1 at the edge of the cell (including the edge of the pattern)');
    pause(2);

%  for non-junctional intProt1 (all but 10 pixels at edge)

    normIntProt1NOjunc = double(intPatt) / double(max(intPatt(:)));
    normIntProt1NOjunc(pattEdge == 1 | roundDistIm < 10) = 0;

    figure; imshow(normIntProt1NOjunc, []);
    title('Non-junctional intProt1');
    pause(2);

% Eventually turn the above into boxplots in R, make sure 10 pixels seems
% right tomorrow!
    
    close all;

%% Storage to structure 'lineOut'

    lineOut(i).filename = full(i).name;
    
    lineOut(i).refProt  = refProt;
    lineOut(i).intProt1 = intProt1;
    lineOut(i).intProt2 = intProt2;
    lineOut(i).dna      = dna;
    
    lineOut(i).pattEdge = pattEdge;
    lineOut(i).pattMask = pattMask;
    lineOut(i).cellMask = cellMask;
    lineOut(i).refWater = refWater;
    
    lineOut(i).refProtThresh     = refProtThresh;
    lineOut(i).refProtThreshEdge = refProtThreshEdge;
    lineOut(i).dnaThresh2x       = dnaThresh2x;

    lineOut(i).normDists   = normDists;
    lineOut(i).roundDists  = roundDists;
    lineOut(i).roundDistIm = roundDistIm;

    lineOut(i).refProtPatt  = double(refPatt) / max(double(refPatt(:)));
    lineOut(i).intProt1Patt = double(intPatt) / max(double(intPatt(:)));

    lineOut(i).intProt1cells = intCells;
    lineOut(i).intProt1normDist = normDistIm;

    lineOut(i).intProt1pattEdge = normIntProt1pattEdge;
    lineOut(i).intProt1junc     = normIntProt1junc;
    lineOut(i).intProt1edge     = normIntProt1edge;
    lineOut(i).intProt1NOjunc   = normIntProt1NOjunc;

end % end loop

%% Vectorizing data across all images and cells for smooth and plot analysis

% Initialize empty arrays for each field
pattMaskVec       = [];
refProtPattVec    = [];
intProt1PattVec   = [];

intProt1cellsVec    = [];
intProt1normDistVec = [];
intProt1pattEdgeVec = [];
intProt1juncVec     = [];
intProt1edgeVec     = [];
intProt1NOjuncVec   = [];

% Loop through all entries in lineOut
for i = 1:numel(lineOut)
    pattMaskVec       = [pattMaskVec;     lineOut(i).pattMask(:)];
    refProtPattVec    = [refProtPattVec;  lineOut(i).refProtPatt(:)];
    intProt1PattVec   = [intProt1PattVec; lineOut(i).intProt1Patt(:)];

    intProt1cellsVec    = [intProt1cellsVec;    lineOut(i).intProt1cells(:)];
    intProt1normDistVec = [intProt1normDistVec; lineOut(i).intProt1normDist(:)];
    intProt1pattEdgeVec = [intProt1pattEdgeVec; lineOut(i).intProt1pattEdge(:)];
    intProt1juncVec     = [intProt1juncVec;     lineOut(i).intProt1junc(:)];
    intProt1edgeVec     = [intProt1edgeVec;     lineOut(i).intProt1edge(:)];
    intProt1NOjuncVec   = [intProt1NOjuncVec;   lineOut(i).intProt1NOjunc(:)];
end

%% Pattern only refProt & intProt1

pattMaskOn       = pattMaskVec == 1;

refProtPattVec   = refProtPattVec(pattMaskOn);
intProt1PattVec  = intProt1PattVec(pattMaskOn);

%% Storing vector variables

lineOut(1).refProtPattVec  = refProtPattVec;
lineOut(1).intProt1PattVec = intProt1PattVec;

lineOut(1).intProt1cellsVec    = intProt1cellsVec;
lineOut(1).intProt1normDistVec = intProt1normDistVec;
lineOut(1).intProt1pattEdgeVec = intProt1pattEdgeVec;
lineOut(1).intProt1juncVec     = intProt1juncVec;
lineOut(1).intProt1edgeVec     = intProt1edgeVec;
lineOut(1).intProt1NOjuncVec   = intProt1NOjuncVec;

%% Create variables for later smooth and plot 

% RefProt vs. intProt1 (intensity only)
refVSintPattVec(:, 1)      = refProtPattVec;
refVSintPattVec(:, 2)      = intProt1PattVec;
lineOut(1).refVSintPattVec = refVSintPattVec(~all(refVSintPattVec == 0, 2), :);

% distance from edge vs. intensity (intProt1 only)
distVintIntProt1(:, 1)      = intProt1normDistVec;
distVintIntProt1(:, 2)      = intProt1cellsVec;
lineOut(1).distVintIntProt1 = distVintIntProt1(~all(distVintIntProt1 == 0, 2), :);

end % end function

%% Helper function for clearing polygons to isolate junctions

function imX = remove_polygon(imX)
    % GUI to draw and remove polygon-defined regions from grayscale image

    % Flags stored in figure UserData
    app.done = false;
    app.draw = false;

    % Create the figure
    hFig = figure('Name', 'Polygon Removal Tool', ...
                  'NumberTitle', 'off', ...
                  'MenuBar', 'none', ...
                  'ToolBar', 'none', ...
                  'Units', 'normalized', ...
                  'Position', [0.2 0.2 0.6 0.7], ...
                  'UserData', app, ...
                  'CloseRequestFcn', @onCloseRequest);

    % Display image
    hAx  = axes('Parent', hFig, 'Position', [0.05 0.25 0.9 0.7]);
    hImg = imshow(imX, [], 'Parent', hAx);
    title(hAx, 'Click "Draw Polygon" to remove a region');

    % Buttons
    uicontrol('Style', 'pushbutton', 'String', 'Draw Polygon', ...
              'Units', 'normalized', ...
              'Position', [0.1 0.05 0.3 0.08], ...
              'FontSize', 12, ...
              'Callback', @(~,~) setFlag('draw', true));

    uicontrol('Style', 'pushbutton', 'String', 'Done', ...
              'Units', 'normalized', ...
              'Position', [0.6 0.05 0.3 0.08], ...
              'FontSize', 12, ...
              'Callback', @(~,~) setFlag('done', true));

    % Main loop
    while ishandle(hFig)
        pause(0.1);
        drawnow;

        app = get(hFig, 'UserData');

        if app.done
            break;
        end

        if app.draw
            title(hAx, 'Draw polygon. Double-click to finish.');
            poly = drawpolygon(hAx, 'Color', 'r');

            if ~isempty(poly) && isvalid(poly)
                BW      = poly.createMask();
                imX(BW) = 0;
                delete(poly);
                set(hImg, 'CData', imX);
                title(hAx, 'Region removed.');
            else
                title(hAx, 'No polygon drawn.');
            end

            setFlag('draw', false);
        end
    end

    if ishandle(hFig)
        delete(hFig);
    end

    % ---- Utility to update flags in UserData ----
    function setFlag(flagName, value)
        app            = get(hFig, 'UserData');
        app.(flagName) = value;
        set(hFig, 'UserData', app);
    end

    % ---- Custom close request ----
    function onCloseRequest(~, ~)
        setFlag('done', true);
    end
end