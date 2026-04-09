function [lineOut] = line_analysis(lineLocker)
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

lineOut(length(lineLocker)) = struct('filename', [], ...
                                     'normAvgIntProt', [], 'avgEucDist', [], ...
                                     'probIntProt', [], 'denIntProt', []);

imInterval = 1:512;
distX      = [];
distY      = [];
intProtInt = [];

%%
for i = 1:length(lineLocker)

    % Pull necessary variables from
    dnaCent       = lineLocker(i).dnaCent;
    refProtThresh = lineLocker(i).refProtThresh;
    intProtNorm   = lineLocker(i).intProtNorm;

    if nnz(refProtThresh) == 0  
    disp(['Skipping iteration ', num2str(i), ' due to all-zero refProtThresh']);
    continue;  % Skip this iteration
    end

    %% Section 1: Load previous alignment locker

    refProtMask = refProtThresh > 0;
    intProtMask = intProtNorm   > 0;

    imshowpair(refProtMask, intProtMask); pause(2);

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
    disp(['Processing file ', num2str(i), ': ', lineLocker(i).filename]);
    disp(['Number of edge pixels: ', num2str(numel(edgeRows))]);

    %% Calculates the distance in only x direction

    % Find the minimum and maximum x-coordinates
    minX = min(edgeCols);
    maxX = max(edgeCols);
    
    % The longest horizontal distance is the difference between minX & maxX
    maxHorizonDist = maxX - minX;
    
    %% Determine y-value from nucleus
        
    intProtYdistNucCent = zeros(512);
    refProtYdistNucCent = zeros(512);

    for y = imInterval
        for x = imInterval
            if intProtMask(y, x) > 0
               intProtYdistNucCent(y, x) = abs(y - dnaCent(2));
            end
            if refProtMask(y, x) > 0
               refProtYdistNucCent(y, x) = abs(y - dnaCent(2));
            end
        end
    end
    
    maxYdistNucCent    = max(refProtYdistNucCent(:));
    relYdistNucCent    = intProtYdistNucCent / maxYdistNucCent;
    relYdistNucCentVec = relYdistNucCent(:);
    distY              = [distY; abs(relYdistNucCentVec)];

    %% Find x distance as a percent of cell length from nuclear centroid for all accepted points in intProt
    
    intProtxDistNucCent = zeros(512);
    refProtxDistNucCent = zeros(512);
    
    for y = imInterval
        for x = imInterval
            if intProtMask(y, x) > 0
               intProtxDistNucCent(y, x) = x - dnaCent(1);
            end
            if refProtMask(y, x) > 0
               refProtxDistNucCent(y, x) = x - dnaCent(1);
            end
        end
    end
    
    maxXdistNucCent    = max(abs(refProtxDistNucCent(:)));
    relXdistNucCent    = intProtxDistNucCent / maxXdistNucCent;
    relXdistNucCentVec = relXdistNucCent(:);
 
    minNegVal = abs(min(relXdistNucCentVec(relXdistNucCentVec < 0))); % Find the absolute value of the minimum negative value
    relXdistNucCent(relXdistNucCent < 0) = relXdistNucCent(relXdistNucCent < 0) / minNegVal;

    for n = 1:length(relXdistNucCent)
        if relXdistNucCentVec(n) < 0
           relXdistNucCentVec(n) = relXdistNucCentVec(n) / minNegVal; 
        end
    end

    distX = [distX; relXdistNucCentVec];

    %% Find euclidian distance as a percent of cell length from nuclear centroid

    eucDistNucCent = zeros(512);

    for y = imInterval
        for x = imInterval
            if intProtMask(y, x) > 0
               eucDistNucCent(y, x) = sqrt((x - dnaCent(1))^2 + (y - dnaCent(2))^2);
               if x < dnaCent(1)
                  eucDistNucCent(y, x) = - eucDistNucCent(y, x);
               end
            end
        end
    end
    
    maxEucDistNucCent = max(abs(eucDistNucCent(:)));
    relEucDistNucCent = eucDistNucCent / maxEucDistNucCent;

    %%

    intProtInt                   = [intProtInt; intProtNorm(:)];
    lineOut(i).filename          = lineLocker(i).filename;
    lineOut(i).intProtnorm       = intProtNorm;
    lineOut(i).relXdistNucCent   = relXdistNucCent;
    lineOut(i).relYdistNucCent   = relYdistNucCent;
    lineOut(i).relEucDistNucCent = relEucDistNucCent;
    lineOut(i).intProtInt        = intProtInt;
    lineOut(i).distX             = distX;
    lineOut(i).distY             = distY;

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
for i = 1:length(lineOut)
    % Get the current cell's data
    currRelEucDistNucCent = lineOut(i).relEucDistNucCent;
    % Append the data to the array
    allEucDist = [allEucDist; currRelEucDistNucCent(:)];  % Convert to column vector if needed
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

lineOut(1).normAvgIntProt = normAvgIntProt;
lineOut(1).avgEucDist     = avgEucDist;

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
   numStacks         = length(lineLocker); % the number of cells used to define normAvgIntProt

   RandProtstack = zeros(size(normAvgIntProt, 1), size(normAvgIntProt, 2), numStacks);
    
   for n = 1:numStacks
       RandProtstack(:, :, n) = rand_protein_sim(normAvgIntProt, 5000, 2);
   end
    
   AvgRandProt = mean(RandProtstack, 3);
    
   normAvgIntProt    = AvgRandProt;
   normAvgIntProtvec = normAvgIntProt(:);
    
   proteinName = 'Random Protein';
end

imagesc(xbinEdges(1:end-1), ybinEdges(1:end-1), normAvgIntProt');
set(gca, 'YDir', 'normal', 'FontName', 'Arial');
colormap(customHot);
axis off;
%colorbar;
%xlabel('Relative X Distance from Nuclear Centroid', 'FontName', 'Arial', 'FontSize', 7);
%ylabel('Relative Y Distance from Nuclear Centroid', 'FontName', 'Arial', 'FontSize', 7);
%title(['Average Normalized ', proteinName, ' Intensity for Average Cell'], 'FontName', 'Arial', 'FontSize', 7);

%%

distVint(:, 1) = avgEucDistVec;
distVint(:, 2) = normAvgIntProtvec;

distVint(distVint(:, 2) == 0, :) = []; % remove rows where intensity is 0

titleStg       = ['Normalized ', proteinName, ' Line Pattern Distance vs Intensity'];
xAxisStg       =  'Normalized Euclidean Distance from Nuclear Midplane';
yAxisStg       = ['Normalized ', proteinName, ' Intensity'];
lineXlim       = [-1 1]; Xtick = -1:0.5:1;
lineYlim       = [ 0 1]; Ytick =  0:0.5:1;

[outDistVint] = smooth_and_plot(distVint, 0.5, titleStg, xAxisStg, yAxisStg, lineXlim, lineYlim, Xtick, Ytick);

lineOut(1).distVint    = distVint;
lineOut(1).outDistVint = outDistVint;

%% 

[sortXdist, sortXindices] = sort(xbinEdges);
sortNormIntProt           = normAvgIntProt(sortXindices);

% Preallocate arrays
sumIndex        = zeros(1, numel(sortXdist));
sumIndexIntProt = zeros(1, numel(sortXdist));

% Total sum of intProt intensities
sumIntProt = sum(sortNormIntProt(:));

% Loop through unique distances
for i = 1:numel(sortXdist)
    x = sortXdist(i);
    
    % Find the last index of x in sortXdist
    lastIdx     = find(sortXdist == x, 1, "last");
    sumIndex(i) = lastIdx;
    
    % Sum up intensities up to this index
    sumIndexIntProt(i) = sum(sortNormIntProt(1:lastIdx));
    ratioSumIntProt(i) = sumIndexIntProt(i) / sumIntProt;

end

% Find the index of the closest value
[~, idx25]  = min(abs(ratioSumIntProt - 0.25));
[~, idx50]  = min(abs(ratioSumIntProt - 0.5));
[~, idx75]  = min(abs(ratioSumIntProt - 0.75));
[~, idx100] = min(abs(ratioSumIntProt - 1));

% Assign probability values based on xbin positions, not indices to intProt

probIntProt = zeros(numXbins, numYbins); % initialize

probIntProt(1:numXbins, xbinEdges >= -1              & xbinEdges <= sortXdist(idx25))  = 0.2;
probIntProt(1:numXbins, xbinEdges > sortXdist(idx25) & xbinEdges <= sortXdist(idx50))  = 0.4;
probIntProt(1:numXbins, xbinEdges > sortXdist(idx50) & xbinEdges <= sortXdist(idx75))  = 0.6;
probIntProt(1:numXbins, xbinEdges > sortXdist(idx75) & xbinEdges <= sortXdist(idx100)) = 0.8;

% Find probIntProt and plot for just cell loc

rotateIntProt = flip(normAvgIntProt', 1);
comboIntProt  = zeros(numXbins, numYbins);
for x = 1:numXbins
    for y = 1:numYbins
        if probIntProt(x, y) && rotateIntProt(x, y) > 0
           comboIntProt(x, y) = probIntProt(x, y);
        end
    end
end

lineOut(1).probIntProt = comboIntProt;

figure; % Initialize new figure window

% Visualize the result
imagesc(comboIntProt);
%title(['Probability Map of ', proteinName], 'FontName', 'Arial', 'FontSize', 10);
colormap(customParula); caxis([0, 1]) % colorbar;
axis off;
hold on;
hold off;

%% Plot density diagram per x bin for IntProt

% Compute bin area
xbinDist = abs(xbinEdges(1) - xbinEdges(2));
ybinDist = abs(ybinEdges(1) - ybinEdges(2));
binArea  = xbinDist * ybinDist;

% Preallocate arrays
binDen     = zeros(size(normAvgIntProt));  % Density per bin
avgDenPerX = zeros(length(normAvgIntProt), 1);  % Avg per x-column

% Preallocate cumulative density array
cumulativeDen = zeros(size(normAvgIntProt));  % Same size as your data

% Compute bin density
for ix = 1:size(normAvgIntProt, 1)
    for iy = 1:size(normAvgIntProt, 2)
        if normAvgIntProt(ix, iy) ~= 0  % Only compute for nonzero bins
           binDen(ix, iy) = normAvgIntProt(ix, iy) / binArea;
        end
    end
end

% Compute average density per x-column (excluding zero-density bins)
for ix = 1:size(binDen, 1)
    nonzeroBins    = binDen(ix, binDen(ix, :) > 0); % Select nonzero bins
    avgDenPerX(ix) = mean(nonzeroBins); % Compute average
end

DenPerX = zeros(size(binDen));
for ix = 1:size(binDen, 1)
    for iy = 1:size(binDen, 2)
        if binDen(ix, iy) > 0
           DenPerX(ix, iy) = avgDenPerX(ix);
        end
    end
end

DenPerX = DenPerX / max(DenPerX(:)); % ormalize

lineOut(1).denIntProt = DenPerX;

% Plot density heatmap
figure;
imagesc(xbinEdges(1:end-1), ybinEdges(1:end-1), DenPerX');  
set(gca, 'YDir', 'normal', 'FontName', 'Arial', 'FontSize', 7);
%title(['Density Map of ', proteinName], 'FontName', 'Arial', 'FontSize', 7);
colormap(customHot); % colorbar;
axis off;