function lineOut = line_confluency_junction_mapping(folder)
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
lineOut(length(full)) = struct(); % NEW

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

   %refProt(refProt < 4500) = 0;
   %intProt1(intProt1 < 5500) = 0;

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

    pattMask1     = refProt > 4500;
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

%% ============================
%% NEW: MANUAL JUNCTION DRAWING
%% ============================

disp('Draw junctions (double click to finish each line, press Done when finished)');

junctionLineMask = draw_junctions_GUI(refProt); % NEW FUNCTION

% Show overlay for validation
figure;
imshow(refProt, []);
hold on;
visboundaries(junctionLineMask, 'Color', 'g');
title('Junction lines overlaid on refProt');
pause(3); % <<< requested 3 sec pause
close;

%% ============================
%% NEW: EXPAND TO JUNCTION REGION (±x pixels)
%% ============================

junctionWidth = 5; % <<< adjustable parameter 'x' as needed

junctionMask = bwdist(junctionLineMask) <= junctionWidth;
junctionMask = junctionMask & pattMask; % NEW (VERY IMPORTANT)

figure;
imshowpair(junctionMask, refProt);
title('Expanded Junction Region (green) over refProt');
pause(2);
close;

%% ============================
%% KEEP: PATTERN MASKING
%% ============================

refPatt = refProt;
refPatt(pattMask == 0) = 0;

intPatt = intProt1;
intPatt(pattMask == 0) = 0;

intPatt2 = intProt2; % NEW
intPatt2(pattMask == 0) = 0; % NEW

%% ============================
%% NEW: NORMALIZATION (PER IMAGE)
%% ============================

%% BACKGROUND SUBTRACTION (NEW)

backThresh = 5000; % or input parameter

intPatt_bs = double(intPatt) - backThresh;
intPatt_bs(intPatt_bs < 0) = 0;

refPatt_bs = double(refPatt) - backThresh;
refPatt_bs(refPatt_bs < 0) = 0;

intPatt2_bs = double(intPatt2) - backThresh;
intPatt2_bs(intPatt2_bs < 0) = 0;

%% NORMALIZATION (AFTER BACKGROUND SUBTRACTION)

int1Max = max(intPatt_bs(:));  if int1Max == 0, int1Max = 1; end
refMax  = max(refPatt_bs(:));  if refMax == 0,  refMax = 1;  end
int2Max = max(intPatt2_bs(:)); if int2Max == 0, int2Max = 1; end

int1Norm = intPatt_bs  / int1Max;
refNorm  = refPatt_bs  / refMax;
int2Norm = intPatt2_bs / int2Max;

%% ============================
%% NEW: DEFINE NON-JUNCTION REGION
%% ============================

nonJunctionMask = pattMask & ~junctionMask;

% OPTIONAL: remove pattern edge (same logic as before)
nonJunctionMask(pattEdge == 1) = 0;
junctionMask(pattEdge == 1)    = 0;

%% ============================
%% NEW: EXTRACT INTENSITIES
%% ============================

junctionVals.intProt1 = int1Norm(junctionMask);
junctionVals.refProt  = refNorm(junctionMask);
junctionVals.intProt2 = int2Norm(junctionMask);

nonJunctionVals.intProt1 = int1Norm(nonJunctionMask);
nonJunctionVals.refProt  = refNorm(nonJunctionMask);
nonJunctionVals.intProt2 = int2Norm(nonJunctionMask);

%% Visualization (sanity check)

figure;
subplot(1,2,1);
imshow(junctionMask, []);
title('Junction Region');

subplot(1,2,2);
imshow(nonJunctionMask, []);
title('Non-Junction Region');

pause(2);
close;

%% ============================
%% STORE (UPDATED)
%% ============================

lineOut(i).junctionMask    = junctionMask;
lineOut(i).nonJunctionMask = nonJunctionMask;

lineOut(i).junction_intProt1    = junctionVals.intProt1;
lineOut(i).junction_refProt     = junctionVals.refProt;
lineOut(i).junction_intProt2    = junctionVals.intProt2;

lineOut(i).nonJunction_intProt1 = nonJunctionVals.intProt1;
lineOut(i).nonJunction_refProt  = nonJunctionVals.refProt;
lineOut(i).nonJunction_intProt2 = nonJunctionVals.intProt2;

end % end loop

%% Vectorizing data across all images and cells for smooth and plot analysis

%% ============================
%% NEW: VECTORIZATION (ALL CHANNELS)
%% ============================

% Initialize
junctionVec_intProt1    = [];
junctionVec_refProt     = [];
junctionVec_intProt2    = [];

nonJunctionVec_intProt1 = [];
nonJunctionVec_refProt  = [];
nonJunctionVec_intProt2 = [];

% Loop through all images
for i = 1:numel(lineOut)

    % --- JUNCTION ---
    junctionVec_intProt1 = [junctionVec_intProt1; lineOut(i).junction_intProt1];
    junctionVec_refProt  = [junctionVec_refProt;  lineOut(i).junction_refProt];
    junctionVec_intProt2 = [junctionVec_intProt2; lineOut(i).junction_intProt2];

    % --- NON-JUNCTION ---
    nonJunctionVec_intProt1 = [nonJunctionVec_intProt1; lineOut(i).nonJunction_intProt1];
    nonJunctionVec_refProt  = [nonJunctionVec_refProt;  lineOut(i).nonJunction_refProt];
    nonJunctionVec_intProt2 = [nonJunctionVec_intProt2; lineOut(i).nonJunction_intProt2];

end

%% ============================
%% STORE (like your previous style)
%% ============================

lineOut(1).junctionVec_intProt1    = junctionVec_intProt1;
lineOut(1).junctionVec_refProt     = junctionVec_refProt;
lineOut(1).junctionVec_intProt2    = junctionVec_intProt2;

lineOut(1).nonJunctionVec_intProt1 = nonJunctionVec_intProt1;
lineOut(1).nonJunctionVec_refProt  = nonJunctionVec_refProt;
lineOut(1).nonJunctionVec_intProt2 = nonJunctionVec_intProt2;

end % end function

%% Helper function for clearing polygons to isolate junctions

function junctionMask = draw_junctions_GUI(im)

    junctionMask = false(size(im));

    % Flags
    app.done = false;
    app.draw = false;

    % Figure
    hFig = figure('Name', 'Draw Junctions', ...
                  'NumberTitle', 'off', ...
                  'MenuBar', 'none', ...
                  'ToolBar', 'none', ...
                  'Units', 'normalized', ...
                  'Position', [0.2 0.2 0.6 0.7], ...
                  'UserData', app);

    hAx = axes('Parent', hFig);
    imshow(im, [], 'Parent', hAx);
    title('Draw junctions (double click to finish each). Click Done when finished.');

    % Buttons
    uicontrol('Style', 'pushbutton', 'String', 'Draw Junction', ...
              'Units', 'normalized', ...
              'Position', [0.1 0.05 0.3 0.08], ...
              'FontSize', 12, ...
              'Callback', @(~,~) setFlag('draw', true));

    uicontrol('Style', 'pushbutton', 'String', 'Done', ...
              'Units', 'normalized', ...
              'Position', [0.6 0.05 0.3 0.08], ...
              'FontSize', 12, ...
              'Callback', @(~,~) setFlag('done', true));

    while ishandle(hFig)
        pause(0.1);
        drawnow;

        app = get(hFig, 'UserData');

        if app.done
            break;
        end

        if app.draw
            h = drawpolyline(hAx, 'Color', 'g');

            if ~isempty(h) && isvalid(h)
                BW = createMask(h);
                junctionMask = junctionMask | BW;
                delete(h);

                % Update display
                imshow(im, [], 'Parent', hAx);
                hold(hAx, 'on');
                visboundaries(junctionMask, 'Color', 'g', 'LineWidth', 0.5);
                hold(hAx, 'off');
            end

            setFlag('draw', false);
        end
    end

    if ishandle(hFig)
        delete(hFig);
    end

    function setFlag(flagName, value)
        app = get(hFig, 'UserData');
        app.(flagName) = value;
        set(hFig, 'UserData', app);
    end
end