function [bowOut] = bow_analysis(bowLocker)
%line analysis:


%% Section 0: Alter color map so that 0 is black

customParula       = parula(256); % Choose an existing colormap (e.g., parula, jet, viridis) and copy it to cmap
customParula(1, :) = [0, 0, 0];   % Set the first entry (corresponding to the lowest value) to black

numColors    = 256;                                                              % Define number of color levels
origColormap = hot(numColors);                                                   % Get the default 'hot' colormap
cutoff       = 0.9;                                                              % Define a cutoff threshold (e.g., 80% of the colormap)
newMaxIdx    = round(cutoff * numColors); 
customHot    = origColormap(1:newMaxIdx, :);                                     % Stretch and cut off colormap
customHot    = [customHot; repmat(customHot(end, :), numColors - newMaxIdx, 1)]; % Extend the last color to fill up to 1

%% Section 1:

bowOut(length(bowLocker)) = struct('filename', [], ...
                                     'normAvgIntProt', [], 'avgEucDist', [], ...
                                     'probIntProt', [], 'denIntProt', []);

imInterval = 1:512;
distX      = [];
distY      = [];
intProtInt = [];

%%
for i = 1:length(bowLocker)

    % Pull necessary variables from
    refCent       = bowLocker(i).refCent;
    refProtThresh = bowLocker(i).refProtThresh;
    intProtNorm   = bowLocker(i).intProtNorm;

    if nnz(refProtThresh) == 0  
    disp(['Skipping iteration ', num2str(i), ' due to all-zero refProtThresh']);
    continue;  % Skip this iteration
    end

    %% Section 1: Load previous alignment locker

    refProtMask = refProtThresh > 0;
    intProtMask = intProtNorm   > 0;

    borderMask                  = zeros(512);
    borderMask(1, imInterval)   = 1;
    borderMask(512, imInterval) = 1;
    borderMask(imInterval, 1)   = 1;
    borderMask(imInterval, 512) = 1;

    [borderRows, borderCols] = find(refProtMask & borderMask == 1);

    edgePixels           = edge(refProtMask, 'sobel'); % Detect edges (can use 'canny' or other methods)
    [edgeRows, edgeCols] = find(edgePixels);         % Get edge pixel coordinates
    
    % Combine the border coordinates with the existing edge coordinates
    edgeRows = [edgeRows; borderRows]; % Append border rows
    edgeCols = [edgeCols; borderCols]; % Append border columns

    %% Calculates euclidean distance-based length of the cell vs the next section that only calculates x distance
  
    % Compute pairwise distances of edge pixels
    dists = pdist([edgeRows, edgeCols]); % Pairwise distances (1D condensed form)
    
    % Convert to squareform distance matrix
    distMat = squareform(dists);
    
    % Find the indices of the points with the longest distance
    [maxDist, linearIdx] = max(distMat(:));                   % Max distance and linear index
    [maxIdx1, maxIdx2]   = ind2sub(size(distMat), linearIdx); % Convert to subscripts
    
    % Plot the results
    imshow(refProtMask); 
    hold on;
    plot(edgeCols, edgeRows, 'r.', 'MarkerSize', 10); % Edge pixels in red
    plot([edgeCols(maxIdx1), edgeCols(maxIdx2)], [edgeRows(maxIdx1), edgeRows(maxIdx2)], 'go-', 'MarkerSize', 12); % Longest distance as green line
    hold off;
    pause(0.5)
    
    % Display the longest distance
    disp(['Longest distance: ', num2str(maxDist)]);
    disp(['Processing file ', num2str(i), ': ', bowLocker(i).filename]);
    disp(['Number of edge pixels: ', num2str(numel(edgeRows))]);

    %% Calculates the distance in only x direction

    % Find the minimum and maximum x-coordinates
    minX = min(edgeCols);
    maxX = max(edgeCols);
    
    % The longest horizontal distance is the difference between minX & maxX
    maxHorizonDist = maxX - minX;
    
    %% Determine y-value from cell centroid
        
intProtYdistCent = zeros(512);
refProtYdistCent = zeros(512);

for y = imInterval
    for x = imInterval
        if intProtMask(y, x) > 0
           intProtYdistCent(y, x) = refCent(2) - y;
        end
        if refProtMask(y, x) > 0
           refProtYdistCent(y, x) = refCent(2) - y;
        end
    end
end

% Normalize with respect to max of reference distances
maxYdistRefCent    = max(abs(refProtYdistCent(:)));
relYdistRefCent    = intProtYdistCent / maxYdistRefCent;

% Flatten to vector
relYdistRefCentVec = relYdistRefCent(:);

% Normalize positives and negatives independently to [0,1] and [0,-1]
maxPos = max(relYdistRefCentVec(relYdistRefCentVec > 0));
minNeg = min(relYdistRefCentVec(relYdistRefCentVec < 0));

% Avoid division by zero
if maxPos ~= 0
    relYdistRefCentVec(relYdistRefCentVec > 0) = ...
        relYdistRefCentVec(relYdistRefCentVec > 0) / maxPos;
end
if minNeg ~= 0
    relYdistRefCentVec(relYdistRefCentVec < 0) = ...
        relYdistRefCentVec(relYdistRefCentVec < 0) / abs(minNeg);
end

% Final values range between -1 and 1
distY = [distY; relYdistRefCentVec];


    %% Find x distance as a percent of cell length from cell centroid for all accepted points in intProt
    
    intProtxDistCent = zeros(512);
    refProtxDistCent = zeros(512);
    
    for y = imInterval
        for x = imInterval
            if intProtMask(y, x) > 0
               intProtxDistCent(y, x) = abs(x - refCent(1));
            end
            if refProtMask(y, x) > 0
               refProtxDistCent(y, x) = abs(x - refCent(1));
            end
        end
    end
    
    maxXdistRefCent    = max(abs(refProtxDistCent(:)));
    relXdistRefCent    = intProtxDistCent / maxXdistRefCent;
    relXdistRefCentVec = relXdistRefCent(:);
 
    distX = [distX; relXdistRefCentVec];

    %% Find euclidian distance as a percent of cell length from cell centroid

    eucDistRefCent = zeros(512);

    for y = imInterval
        for x = imInterval
            if intProtMask(y, x) > 0
               eucDistRefCent(y, x) = sqrt((x - refCent(1))^2 + (y - refCent(2))^2);
               if y < refCent(2)
                  eucDistRefCent(y, x) = - eucDistRefCent(y, x);
               end
            end
        end
    end
    
    maxEucDistRefCent = max(abs(eucDistRefCent(:)));
    relEucDistRefCent = eucDistRefCent / maxEucDistRefCent;

    %%

    intProtInt                  = [intProtInt; intProtNorm(:)];
    bowOut(i).filename          = bowLocker(i).filename;
    bowOut(i).intProtnorm       = intProtNorm;
    bowOut(i).relXdistRefCent   = relXdistRefCent;
    bowOut(i).relYdistRefCent   = relYdistRefCent;
    bowOut(i).relEucDistRefCent = relEucDistRefCent;
    bowOut(i).intProtInt        = intProtInt;
    bowOut(i).distX             = distX;
    bowOut(i).distY             = distY;

end
%% 

% Define bin edges for distances
xbinEdges = linspace(min(distX), max(distX), 51);  % bin edges for x distances
ybinEdges = linspace(min(distY), max(distY), 51);  % bin edges for y distances
numXbins  = length(xbinEdges) - 1;
numYbins  = length(ybinEdges) - 1;

% Initialize arrays to store average intensities and distances for each bin
avgIntProtInt = zeros(numXbins, numYbins);
avgEucDist    = zeros(numXbins, numYbins);

% Now we need to build the array for `currEucDist` similar to how
% we did for `distances_x` and `distances_y`
allEucDist = [];

% Loop over all cells to create the `currEucDist` equivalent for each
for i = 1:length(bowOut)
    % Get the current cell's data
    currRelEucDistRefCent = bowOut(i).relEucDistRefCent;
    % Append the data to the array
    allEucDist = [allEucDist; currRelEucDistRefCent(:)];  % Convert to column vector if needed
end

% Loop through each bin and calculate the average intensity and distance
for ix = 1:numXbins
    for iy = 1:numYbins
        % Get indices of points within the current bin range
        idxX = distX >= xbinEdges(ix) & distX < xbinEdges(ix+1);
        idxY = distY >= ybinEdges(iy) & distY < ybinEdges(iy+1);
        
        % Find common indices for both x and y bin ranges
        validIdx = idxX & idxY;
        
        if any(validIdx)
           % Calculate the average intensity and distance for these points
           avgIntProtInt(ix, iy) = mean(intProtInt(validIdx));
           avgEucDist(ix, iy)    = mean(allEucDist(validIdx));  % Use `allEucDist`
        end
    end
end

normAvgIntProt = avgIntProtInt / max(avgIntProtInt(:));

bowOut(1).normAvgIntProt = normAvgIntProt;
bowOut(1).avgEucDist     = avgEucDist;

%% Section 1: Pull necessary variables from lineOut structure 

normAvgIntProtvec = normAvgIntProt(:);
avgEucDistVec     = avgEucDist(:);

%%

% Prompt user to choose between Random Protein or Protein of Interest
proteinChoice = questdlg(...
    'Would you like to analyze a Random Protein or a Protein of Interest?', ...
    'Protein Analysis Option', ...
    'Random Protein', 'Protein of Interest', 'Random Protein');

% Handle case where user closes dialog
if isempty(proteinChoice)
   error('No option selected. Exiting script.');
end

%%
% Plot the average IntProt intensity as a heatmap

if strcmp(proteinChoice, 'Protein of Interest')
   proteinName = inputdlg(...
       'What is the name of your protein of interest?', ...
       'Protein Name Input', ...
       [1 50]);  % One-line input, 50-character width

   % Handle case where user cancels or leaves input empty
   if isempty(proteinName) || isempty(proteinName{1})
      error('No protein name entered. Exiting script.');
   end

   % Extract the actual string
   proteinName = proteinName{1};

else
   
   normAvgIntProtvec = []; 
   numStacks         = length(bowLocker); % the number of cells used to define normAvgIntProt

   RandProtstack = zeros(size(normAvgIntProt, 1), size(normAvgIntProt, 2), numStacks);
    
   for n = 1:numStacks
       RandProtstack(:, :, n) = rand_protein_sim(normAvgIntProt, 50, 2);
   end
    
   AvgRandProt = mean(RandProtstack, 3);
    
   normAvgIntProt    = AvgRandProt;
   normAvgIntProtvec = normAvgIntProt(:);
    
   proteinName = 'Random Protein';
end

imagesc(xbinEdges(1:end-1), ybinEdges(1:end-1), normAvgIntProt');
set(gca, 'YDir', 'normal', 'FontName', 'Arial');
colormap(customHot);
axis tight;
%colorbar;
xlabel('Relative X Distance from Cell Centroid', 'FontName', 'Arial', 'FontSize', 7);
ylabel('Relative Y Distance from Cell Centroid', 'FontName', 'Arial', 'FontSize', 7);
title(['Average Normalized ', proteinName, ' Intensity for Average Cell'], 'FontName', 'Arial', 'FontSize', 7);

%%

distVint(:, 1) = avgEucDistVec;
distVint(:, 2) = normAvgIntProtvec;

distVint(distVint(:, 2) == 0, :) = []; % remove rows where intensity is 0

titleStg       = ['Normalized ', proteinName, ' Line Pattern Distance vs Intensity'];
xAxisStg       =  'Normalized Euclidean Distance from Cell Midplane';
yAxisStg       = ['Normalized ', proteinName, ' Intensity'];
lineXlim       = [-1 1]; Xtick = -1:0.5:1;
lineYlim       = [ 0 1]; Ytick =  0:0.5:1;

[outDistVint] = smooth_and_plot(distVint, 0.5, titleStg, xAxisStg, yAxisStg, lineXlim, lineYlim, Xtick, Ytick);

bowOut(1).distVint    = distVint;
bowOut(1).outDistVint = outDistVint;

%% 

[sortYdist, sortYindices] = sort(ybinEdges(1:end-1));
sortNormIntProt           = normAvgIntProt(:, sortYindices);

% Preallocate arrays
sumIndex        = zeros(1, numel(sortYdist));
sumIndexIntProt = zeros(1, numel(sortYdist));
ratioSumIntProt = zeros(1, numel(sortYdist));

% Total sum of protein intensities
sumIntProt = sum(sortNormIntProt(:));

% Loop through sorted Y positions
for i = 1:numel(sortYdist)
    y = sortYdist(i);

    % Find last index with value y
    lastIdx     = find(sortYdist == y, 1, "last");
    sumIndex(i) = lastIdx;

    % Sum intensities up to that column
    sumIndexIntProt(i) = sum(sortNormIntProt(:, 1:lastIdx), 'all');
    ratioSumIntProt(i) = sumIndexIntProt(i) / sumIntProt;
end

% Find index for quartiles
[~, idx25]  = min(abs(ratioSumIntProt - 0.25));
[~, idx50]  = min(abs(ratioSumIntProt - 0.5));
[~, idx75]  = min(abs(ratioSumIntProt - 0.75));
[~, idx100] = min(abs(ratioSumIntProt - 1));

% Create and assign probability values based on ybinEdges (along rows)
probIntProt = zeros(numXbins, numYbins);

probIntProt(ybinEdges >= -1              & ybinEdges <= sortYdist(idx25),  1:numXbins) = 0.2;
probIntProt(ybinEdges > sortYdist(idx25) & ybinEdges <= sortYdist(idx50),  1:numXbins) = 0.4;
probIntProt(ybinEdges > sortYdist(idx50) & ybinEdges <= sortYdist(idx75),  1:numXbins) = 0.6;
probIntProt(ybinEdges > sortYdist(idx75) & ybinEdges <= sortYdist(idx100), 1:numXbins) = 0.8;

% Transpose and flip intensity map to match plotting orientation
rotateIntProt = flip(normAvgIntProt', 1);

% Multiply only where both intensity and probability are non-zero
comboIntProt = zeros(numXbins, numYbins);
for x = 1:numXbins
    for y = 1:numYbins
        if probIntProt(x, y) && rotateIntProt(x, y) > 0
            comboIntProt(x, y) = probIntProt(x, y);
        end
    end
end

% Store output
bowOut(1).probIntProt = comboIntProt;

% Plot result
figure;
imagesc(comboIntProt);
title(['Probability Map of ', proteinName], 'FontName', 'Arial', 'FontSize', 10);
colormap(customParula); caxis([0, 1]);
axis off;

%% Plot density diagram per x bin for IntProt

% Compute bin area
xbinDist = abs(xbinEdges(1) - xbinEdges(2));
ybinDist = abs(ybinEdges(1) - ybinEdges(2));
binArea  = xbinDist * ybinDist;

% Preallocate arrays
binDen     = zeros(size(normAvgIntProt));           % Density per bin
avgDenPerY = zeros(size(normAvgIntProt, 2), 1);      % Avg per y-row

% Compute bin density (intensity per unit area)
binDen(normAvgIntProt ~= 0) = normAvgIntProt(normAvgIntProt ~= 0) / binArea;

% Compute average density per y-row, excluding zero bins
for iy = 1:size(binDen, 2)
    nonzeroBins = binDen(binDen(:, iy) > 0, iy);
    if ~isempty(nonzeroBins)
        avgDenPerY(iy) = mean(nonzeroBins);
    end
end

% Fill DenPerY with average value per y-row, where original data is nonzero
DenPerY = zeros(size(binDen));
for iy = 1:size(binDen, 2)
    for ix = 1:size(binDen, 1)
        if binDen(ix, iy) > 0
           DenPerY(ix, iy) = avgDenPerY(iy);
        end
    end
end

% Normalize the density map
DenPerY = DenPerY / max(DenPerY(:));

% Store result
bowOut(1).denIntProt = DenPerY;

% Plot density heatmap
figure;
imagesc(xbinEdges(1:end-1), ybinEdges(1:end-1), DenPerY');  
set(gca, 'YDir', 'normal', 'FontName', 'Arial', 'FontSize', 7);
title(['Density Map of ', proteinName], 'FontName', 'Arial', 'FontSize', 7);
colormap(customHot); % colorbar;
axis off;