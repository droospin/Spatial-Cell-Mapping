function mouseLocker = mouseModelAnalysis(folder)

choice = questdlg('Are you modeling a stable vessel or a migratory front?', ...
                  'Stable Vessel or Migratory Front?', 'Stable Vessel','Migratory Front', 'Migratory Front');

%% Section 1: Get list of all TIF files in the folder, including subfolders

% Get list of all TIF files in the folder, including subfolders
imFiles = fullfile(folder, '**', '*.nd2');
full    = dir(imFiles);

%% Main Loop

for i = 1:length(full)

    fprintf('%d out of %d\n', i, length(full)); % display image number compared to total
    disp(full(i).name); % display image title
    
    %% Image reading
    source = fullfile(full(i).folder, full(i).name);
    
    vesselData = read_nd2(source);

    dna     = vesselData.dna;
    refProt = vesselData.refProt;
    intProt = vesselData.intProt;

    %%

    if strcmp(choice, 'Stable Vessel')
    
       
    title_stg = 'Displaying raw channels. Draw a polygon to KEEP the region.';
    h         = rgbToggleGUI(refProt, intProt, dna, title_stg);
    
    % Get polygon position
    poly_pos = h.Position;
    
    % Create a binary mask from polygon
    polyMask = roipoly(refProt, poly_pos(:,1), poly_pos(:,2));
    
    % Apply mask to each channel
    refProtThresh1            = double(refProt);
    refProtThresh1(~polyMask) = 0;
    intProtThresh1            = double(intProt);
    intProtThresh1(~polyMask) = 0;
    dnaThresh                 = double(dna);
    dnaThresh(~polyMask)      = 0;
    
    % Display result
    figure;
    imshowpair(refProtThresh1, intProtThresh1, 'montage');
    title('Polygon region kept.');
    pause(2);
    
    % Ask if user is happy
    polyChoice = questdlg('Are you happy with this selection?', ...
                      'Confirm Region', ...
                      'Yes', 'No - Redraw', 'Yes');
    
    switch polyChoice
        case 'Yes'
            % Continue with your code
            disp('User confirmed selection.');
            
        case 'No - Redraw'
            disp('User wants to redraw. Restarting selection...');
            close all;
            % You can wrap the code block in a while-loop to allow redrawing:
            % e.g., set a flag: keepGoing = false; then repeat from the top
    end

    %% Remove non-cellular region
    
    else

    title_stg = 'Displaying raw channels. Draw a polygon to REMOVE the region (from migratory front).';
    h         = rgbToggleGUI(refProt, intProt, dna, title_stg);
    
    % Get polygon position
    poly_pos = h.Position;
    
    % Create a binary mask from polygon
    polyMask = roipoly(refProt, poly_pos(:,1), poly_pos(:,2));
    
    % Zero out pixels inside polygon
    refProtThresh1           = double(refProt);
    refProtThresh1(polyMask) = 0;
    intProtThresh1           = double(intProt);
    intProtThresh1(polyMask) = 0;
    dnaThresh                = double(dna);
    dnaThresh(polyMask)      = 0;
    
    figure;
    imshowpair(refProtThresh1, intProtThresh1, 'montage'); 
    title('Polygon region removed.'); pause(2);

    end
    
    %%
    
    refProtThresh2                      = refProtThresh1;
    refProtThresh2(refProtThresh1 < 50) = 0; 
    refProtMask                         = refProtThresh2 > 0;
    imshow(refProtMask); title('Display of crude intensity threshold of refProt'); pause(2);
    
    %%
    
    minSizeThresh      = 100;
    refProtRegionprops = regionprops(refProtMask, 'Area');                                   % Compute region properties (Area) for connected components in the binary mask
    validRegions       = [refProtRegionprops.Area] >= minSizeThresh;  
    refProtMask        = ismember(labelmatrix(bwconncomp(refProtMask)), find(validRegions)); % Create a mask that includes only valid connected components
    
    refProtThresh2(refProtMask == 0) = 0;
    intProtThresh2                   = intProtThresh1;
    intProtThresh2(refProtMask == 0) = 0;
    
    dnaThresh1                   = dnaThresh;
    dnaThresh1(refProtMask == 0) = 0;
    dnaFill                      = imfill(dnaThresh1, 'holes');
    dnaThresh1                   = dnaThresh;
    dnaThresh1(~dnaFill)         = 0;
    
    imshowpair(refProtThresh2, intProtThresh2); title('Display size thresholding of all proteins')
    pause(2); close all;
    
    %% Repeatable polygon-based cell removal (non-endothelial)
    
    intProtEC = intProtThresh2;  % Start with thresholded images
    refProtEC = refProtThresh2;
    
    done = false;
    
    while ~done
           title_stg = 'Draw polygon to REMOVE region. Double-click to finish.';
           p         = rgbToggleGUI(refProtEC, intProtEC, dnaThresh1, title_stg);
    
           % Skip if user deletes the ROI or cancels
           if isempty(p) || isempty(p.Position)
              disp('No polygon drawn. Skipping...');
              break;
           end
    
           % Get polygon coordinates
           poly_pos = p.Position;
        
           % Create binary mask from polygon
           polyMask = roipoly(intProtEC, poly_pos(:,1), poly_pos(:,2));
        
           % Zero out pixels inside polygon
           intProtEC(polyMask) = 0;
           refProtEC(polyMask) = 0;
    
           close;  % Close current figure
    
           % Ask user if they want to draw another
           removeChoice = questdlg('Remove another region?', ...
               'Continue?', ...
               'Yes', 'No', 'No');
    
           if strcmp(removeChoice, 'No')
              done = true;
           end
    end

    close all;
    % Show result
    figure;
    imshowpair(refProtEC, intProtEC, 'montage');
    title('Final result: all marked regions removed');
    close all;
    
    %% Draw line to indicate migratory/angiogenic front
    
    if strcmp(choice, 'Stable Vessel')
       
       disp("No need to draw migratory line because you aren't observing the migratory front.")

       %% Save variables

       mouseLocker(i).filename   = full(i).name;
       mouseLocker(i).vesselType = choice;

       mouseLocker(i).refProt   = refProt;
       mouseLocker(i).refProtEC = refProtEC;

       mouseLocker(i).intProt   = intProt;
       mouseLocker(i).intProtEC = intProtEC;

       mouseLocker(i).dna       = dna;
       mouseLocker(i).dnaThresh = dnaThresh1;

    else

    title_stg = 'Draw line to indicate migratory/angiogenic front.';
    l         = rgbToggleGUI(refProtEC, intProtEC, dnaThresh1, title_stg);
    
    % Get vertex coordinates
    vertices = l.Position;
    
    % Create a blank RGB canvas
    canvasRGB = zeros([size(refProtThresh2), 3], 'uint8');
    
    % Line thickness (in pixels) — adjust this as needed
    lineWidth = 1;
    
    % Loop through each pair of points in the polyline and draw lines
    for n = 1:size(vertices, 1) - 1
        p1 = vertices(n, :);
        p2 = vertices(n + 1, :);
        canvasRGB = insertShape(canvasRGB, 'Line', [p1 p2], ...
            'Color', 'white', 'LineWidth', lineWidth, 'SmoothEdges', true);
    end
    
    % Convert to grayscale and threshold to generate binary mask
    grayLine = rgb2gray(canvasRGB);
    mask     = grayLine > 0;  % or adjust threshold for softer edges if desired
    
    close all;
    
    % Optionally: display result
    imshowpair(refProtThresh2, mask); 
    title('Display of drawn migratory front line'); pause(2);

    
    %% Calculate distance from non-zero points to migratory front for refProt & intProt
       
    % Compute an “all‑pixels” distance map to that line
    distMap       = bwdist(mask);  % Euclidean distance (in pixels) to nearest line pixel
    distMap(mask) = 0;             % exactly zero on the line itself
    
    % Copy the distances just for the positive‑signal pixels
    refProtDist                = zeros(size(refProtEC));   % or size(refProtThresh)
    intProtDist                = zeros(size(intProtEC));
    refProtDist(refProtEC > 0) = distMap(refProtEC > 0);
    intProtDist(intProtEC > 0) = distMap(intProtEC > 0);
    
    % (Optional) visualise as a heat map
    figure;
    imshow(refProtDist, [], 'InitialMagnification', 'fit');
    colormap hot; colorbar;
    title('Distance (px) from migratory/angiogenic front');
    close all;
    
    mouseLocker(i).filename   = full(i).name;
    mouseLocker(i).vesselType = choice;

    mouseLocker(i).refProt     = refProt;
    mouseLocker(i).refProtEC   = refProtEC;
    mouseLocker(i).refProtDist = refProtDist;

    mouseLocker(i).intProt     = intProt;
    mouseLocker(i).intProtEC   = intProtEC;
    mouseLocker(i).intProtDist = intProtDist;

    mouseLocker(i).dna       = dna;
    mouseLocker(i).dnaThresh = dnaThresh1;

    end % choice question dialog
end % main loop

    if strcmp(choice, 'Stable Vessel')

       refProtECvec   = [];
       intProtECvec   = [];

       for i = 1:numel(mouseLocker)
           
           refProtECvec = [refProtECvec; mouseLocker(i).refProtEC(:)];
           intProtECvec = [intProtECvec; mouseLocker(i).intProtEC(:)];
           idx          = refProtECvec > 0;

       end

       mouseLocker(1).refProtVintProt = [refProtECvec(idx), intProtECvec(idx)];

    else
    
    %% Vectorize distances and intensities

    refProtDistVec = [];
    intProtDistVec = [];
    refProtECvec   = [];
    intProtECvec   = [];

    for i = 1:numel(mouseLocker)
    
        refProtDistVec = [refProtDistVec; mouseLocker(i).refProtDist(:)];
        intProtDistVec = [intProtDistVec; mouseLocker(i).intProtDist(:)];
    
        refProtECvec = [refProtECvec; mouseLocker(i).refProtEC(:)];
        intProtECvec = [intProtECvec; mouseLocker(i).intProtEC(:)];
    
    end

    %% Create distance vs. intensity variables

    idx = refProtECvec > 0;
    mouseLocker(1).refProtDistVInt = [refProtDistVec(idx), refProtECvec(idx)];

    idx = intProtECvec > 0;
    mouseLocker(1).intProtDistVInt = [intProtDistVec(idx), intProtECvec(idx)];

    %% Create intProt vs. refProt variables

    mouseLocker(1).refProtVintProt = [refProtECvec(idx), intProtECvec(idx)];

    end

end % function

%% Read image data

function vesselData = read_nd2(nd2file)
    % Read a single .nd2 file and extract protein channels

    imObj = BioformatsImage(nd2file);

    vesselData.filename = nd2file;
    vesselData.xSize    = imObj.height;
    vesselData.ySize    = imObj.width;
    
    vesselData.dna     = getPlane(imObj, 1, 'Cy5', 1);   % dna protein
    vesselData.refProt = getPlane(imObj, 1, 'FITC', 1);  % Reference protein
    vesselData.intProt = getPlane(imObj, 1, 'TRITC', 1); % Interested protein
     
end

%%

function lineORshape = rgbToggleGUI(redChannel, greenChannel, blueChannel, title_stg)

    close all;
    
    % Normalize channels
    redNorm   = mat2gray(redChannel);
    greenNorm = mat2gray(greenChannel);
    blueNorm  = mat2gray(blueChannel);

    % Set up figure size
    figW = size(redNorm, 2) + 40;  % Add padding
    figH = size(redNorm, 1) + 200;
    hFig = figure('Name', 'RGB Channel Toggle', 'NumberTitle', 'off', ...
                  'Position', [100 100 figW figH]);
    
    % Set up axes
    hAx = axes('Parent', hFig, 'Units', 'pixels', ...
               'Position', [20 100 size(redNorm,2) size(redNorm,1)]);  % Leave room at bottom
    
    % Create a static text uicontrol below the image
    hTitleLabel = uicontrol('Style', 'text', ...
                            'Parent', hFig, ...
                            'String', title_stg, ...
                            'Units', 'pixels', ...
                            'Position', [20 60 size(redNorm, 2) 20], ...
                            'HorizontalAlignment', 'center', ...
                            'FontWeight', 'bold');

    % Create checkboxes for each channel
    cbRed = uicontrol('Style', 'checkbox', 'String', 'Red Channel', ...
                      'Value', 1, 'Position', [20 70 120 30], ...
                      'Callback', @(src, event) updateDisplay());
    cbGreen = uicontrol('Style', 'checkbox', 'String', 'Green Channel', ...
                        'Value', 1, 'Position', [20 40 120 30], ...
                        'Callback', @(src, event) updateDisplay());
    cbBlue = uicontrol('Style', 'checkbox', 'String', 'Blue Channel', ...
                       'Value', 1, 'Position', [20 10 120 30], ...
                       'Callback', @(src, event) updateDisplay());

    % Add buttons for drawing
    btnLine = uicontrol('Style', 'pushbutton', 'String', 'Draw Polyline', ...
                        'Position', [160 70 120 30], ...
                        'Callback', @(src, event) drawShape('polyline'));
    btnPoly = uicontrol('Style', 'pushbutton', 'String', 'Draw Polygon', ...
                        'Position', [160 30 120 30], ...
                        'Callback', @(src, event) drawShape('polygon'));

    % Initialize output
    lineORshape = [];

    % Initial display
    updateDisplay();

    % Wait for user to draw and close
    uiwait(hFig);

    % ---------- Nested functions ----------
    function updateDisplay()
        % Get checkbox states
        showRed   = get(cbRed, 'Value');
        showGreen = get(cbGreen, 'Value');
        showBlue  = get(cbBlue, 'Value');

        % Create RGB image
        R               = zeros(size(redNorm)); G = R; B = R;
        if showRed,   R = redNorm;   end
        if showGreen, G = greenNorm; end
        if showBlue,  B = blueNorm;  end

        rgbImage = cat(3, R, G, B);

        % Display with fixed magnification and title
        imshow(rgbImage, 'Parent', hAx, 'InitialMagnification', 100);
        title(hAx, title_stg, 'FontSize', 12, 'Interpreter', 'none');
    end

    function drawShape(type)
        axes(hAx); % Activate axes before drawing
        switch type
            case 'polyline'
                h = drawpolyline('Parent', hAx);
            case 'polygon'
                h = drawpolygon('Parent', hAx);
        end

        if isvalid(h)
            lineORshape = h;
            uiresume(hFig); % Return to calling function after draw
        end
    end
end