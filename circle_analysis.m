function [circleOut] = circle_analysis(circleLocker)
%circle_analysis: 

%% Section 0: Alter color map so that 0 is black

customParula       = parula(256); % Choose an existing colormap (e.g., parula, jet, viridis) and copy it to cmap
customParula(1, :) = [0, 0, 0];   % Set the first entry (corresponding to the lowest value) to black


numColors    = 256;                                                              % Define number of color levels
origColormap = hot(numColors);                                                   % Get the default 'hot' colormap
cutoff       = 0.9;                                                              % Define a cutoff threshold (e.g., 80% of the colormap)
newMaxIdx    = round(cutoff * numColors); 
customHot    = origColormap(1:newMaxIdx, :);                                     % Stretch and cut off colormap
customHot    = [customHot; repmat(customHot(end, :), numColors - newMaxIdx, 1)]; % Extend the last color to fill up to 1

%% Section 1: Develop an all image z-stack 

imInterval   = 1:512;                                 % x- & y-interval for our images (will be used later in code as well)
intProtstack = zeros(512, 512, length(circleLocker)); % initialize intProtstack 

for i = 1:length(circleLocker)
    intProtstack(:, :, i) = circleLocker(i).intProttranslate;
end


%randProtstack = rand_protein_sim(intProtstack, 5000, 2);


circleOut(1).intProtstack  = intProtstack;

%% Section 2: Create 'average cell' based on all the standardized images from max_stack_intProt across the 3rd dimension

% Count valid (nonzero) entries per (x, y) position
nonzeroCounts = sum(intProtstack > 0, 3); % Works for both intProtstack and randProtstack

% Sum nonzero values across slices
sumIntProt  = sum(intProtstack, 3);

% Avoid division by zero
nonzeroCounts(nonzeroCounts == 0) = 1;

% Compute average only for nonzero values
avgCellProt = sumIntProt ./ nonzeroCounts;

%% Make User Decide between random or real protein imaging
% Prompt user to choose between Random Protein or Protein of Interest
proteinChoice = questdlg(...
    'Would you like to analyze a Random Protein or a Protein of Interest?', ...
    'Protein Analysis Option', ...
    'Random Protein', 'Protein of Interest', 'Random Protein'); % This is where you tell MatLab what data to present (real or randomized)

% Handle case where user closes dialog
if isempty(proteinChoice)
   error('No option selected. Exiting script.'); %It exits script if you dont decide between the two paths
end

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
   proteinName = proteinName{1}; % Calls on the name you assigned the protein

   backThresh = 3500;  % subtract background based on image (2750 for actin; 4000 for cav1 for Gen1 Circle Images)

else
%Ensures you are generating a sufficient number of randomized images so
%accurate parallels can be drawn
numStacks     = length(circleLocker); 
randProtstack = zeros(512, 512, numStacks);
%Make a stack of random images based on number of real images:
for n = 1:numStacks
    randProtstack(:, :, n) = rand_protein_sim(avgCellProt, 5000, 2);
end

%Avg all that:
nonzeroCounts = sum(randProtstack > 0, 3);
nonzeroCounts(nonzeroCounts == 0) = 1;

sumRandProt = sum(randProtstack, 3);
avgCellProt = sumRandProt ./ nonzeroCounts; %this is the same as finding a mean of your randomized protein images
%Save it:
circleOut(1).randProtstack = randProtstack;
proteinName = 'Random'; %When you choose to make figures with the randomized proteins

backThresh = 0; 

end

%% Section 3: Baseline thresholding of our average cell

avgCellIntProtThresh = zeros(512, 512); % Initialize average_cell_thresh_intProt

for x = imInterval
    for y = imInterval
        if avgCellProt(x, y) >= backThresh                              % If our average cell is above the background threshold set above
           avgCellIntProtThresh(x, y) = avgCellProt(x, y) - backThresh; % then store that value minus the background threshold into avgCellIntProtThresh
        end
    end
end

imshow(avgCellIntProtThresh, []); % inspect our autoscaled ([]) average cell

%% Section 4: Discount small regions interfering with clean average (areas not part of entire image)

minSizeThresh    = 100;                      % You can change this value as per your requirement. A standard HUVEC cell should be at least 8000, though
binaryAvgIntProt = avgCellIntProtThresh > 0;

intProtregionprops = regionprops(binaryAvgIntProt, 'Area'); % Compute properties (for area and pixel locations) of the segmented objects

for i = 1:numel(intProtregionprops)                % Update areas based on the filtered values.
    if intProtregionprops(i).Area < minSizeThresh  % if Area(i) is smaller than threshold,
       intProtregionprops(i).Area = 0;             % set equal to zero
    end
end

avgCellIntProtMask = ismember(labelmatrix(bwconncomp(binaryAvgIntProt)), find([intProtregionprops.Area] > 0)); % Generate a binary image with filtered objects

avgCell2xThreshIntProt = zeros(512);
for x = imInterval
    for y = imInterval
        if avgCellIntProtMask(x, y) == 1
           avgCell2xThreshIntProt(x, y) = avgCellIntProtThresh(x, y);
        end
    end
end

imshow(avgCell2xThreshIntProt, []);

%% Section 5: Calculating the average cell's distance from the center for all it's points 

binaryAvgIntProt   = avgCell2xThreshIntProt > 0;                % create binary version of our average cell for use in regionprops
intProtregionprops = regionprops(binaryAvgIntProt, 'Centroid'); % calculate the centroid of our binary average cell
intProtCent        = intProtregionprops(1).Centroid;            % grab the most central centroid from the structure 
                                                                % (there might be small features that are picked up on the periphery)

unscaledIntProtdistFromCent = zeros(512, 512); % initialize unscaledIntProtdistFromCent

for x = imInterval
    for y = imInterval
        if binaryAvgIntProt(x, y) > 0
           unscaledIntProtdistFromCent(x, y) = sqrt((x - intProtCent(1))^2 + (y - intProtCent(2))^2); % Calculate the Euclidean distance from (x, y) to the centroid
        end
    end
end

% Visualize
imshow(unscaledIntProtdistFromCent, []);
hold on;
%colorbar;
plot(intProtCent(1), intProtCent(2), 'rx', 'MarkerSize', 10, 'LineWidth', 2); % Ensure you are calculating distances from the proper centroid
hold off;

%% Section 6: Calculating the average cell's distance from the nearest edge for all it's points 

binaryZeroMask          = (binaryAvgIntProt == 0); % Create a binary mask for subsequent use of bwdist
unscaledIntProtdistEdge = bwdist(binaryZeroMask);  % Compute the distance to the nearest zero point (edge point)

% Visualize
imshow(unscaledIntProtdistEdge, []);              % Ensure your calculation worked properly
colorbar;

%% Section 7: Standardizing our previous measurements in section 5 and 6 from pixels to microns

% Define the conversion factor from pixels to microns
pixel2um60x = 0.1872;               

% Scale the distance measurements from Section 6 to microns
scaledIntProtdistFromCent = unscaledIntProtdistFromCent * pixel2um60x; 
scaledIntProtdistFromEdge = unscaledIntProtdistEdge     * pixel2um60x;

% Round up the scaled distances to the nearest integer
roundScaledIntProtdistFromCent = ceil(scaledIntProtdistFromCent); 
roundScaledIntProtDistFromEdge = ceil(scaledIntProtdistFromEdge);

% Determine the maximum distance for center and edge (in microns)
maxDistFromCent = max(scaledIntProtdistFromCent(:)); 
maxDistFromEdge = max(scaledIntProtdistFromEdge(:));

% Find the row and column indices of non-zero elements in the center distance matrix
[rows, cols] = find(binaryAvgIntProt);

% Extract the corresponding distance values for these points
unscaledIntProtDistPtsFromCent = unscaledIntProtdistFromCent(sub2ind(size(unscaledIntProtdistFromCent), rows, cols)); 
unscaledIntProtDistPtsFromEdge = unscaledIntProtdistEdge(sub2ind(size(unscaledIntProtdistEdge), rows, cols));

% Scale the extracted distance points to microns
intProtdistFromCenter = unscaledIntProtDistPtsFromCent * pixel2um60x; 
intProtdistFromEdge   = unscaledIntProtDistPtsFromEdge * pixel2um60x;

% Normalize the scaled distances by dividing by their respective maximums
normIntProtdistFromCent = intProtdistFromCenter / maxDistFromCent; 
normIntProtdistFromEdge = intProtdistFromEdge   / maxDistFromEdge;

% Display the rounded, scaled center distances for verification
imshowpair(roundScaledIntProtdistFromCent, roundScaledIntProtDistFromEdge, 'montage');

%% Section 8: Vectorizing and normalizing the average thresholded cell

intProtInt = avgCell2xThreshIntProt(sub2ind(size(avgCell2xThreshIntProt), rows, cols)); % store all the intensity values of our average cell into a 1D vector
maxInt     = max(avgCell2xThreshIntProt(:));                                            % Get the max fluorescent intensity for the image for later normalization

normIntProtInt = intProtInt / maxInt;             % Normalize fluorescent intensities from our 1D vector
normAvgIntProt = avgCell2xThreshIntProt / maxInt; % Normalize our image of the average cell

% Create heat map for the normalized average intProt
figure;
imagesc(normAvgIntProt);
colormap(customHot)
%colorbar;
axis off;
%title(['Average ', proteinName, ' Cell']);

% Get exactly what MATLAB is displaying:
F      = getframe(gca);      % exact pixels
RGB    = F.cdata;            % exact colors
imwrite(RGB, 'normAvgIntProt.png');  % identical in appearance

%% Section 12:

distVint(:, 1) = normIntProtdistFromEdge;
distVint(:, 2) = normIntProtInt;
titleStg       = [proteinName, ' Distance vs Intensity']; %the protein name you assigned will be given in the title
xAxisStg       = 'Normalized Distance';
yAxisStg       = 'Normalized Intensity';
circleXlim     = [0 1]; Xtick = 0:0.2:1;
circleYlim     = [0 1]; Ytick = 0:0.2:1;

[outDistVint] = smooth_and_plot(distVint, 0.05, titleStg, xAxisStg, yAxisStg, circleXlim, circleYlim, Xtick, Ytick);
circleOut.distVint = distVint;

%% Section 13: IntProt Probability calculations and plotting

roundDist                 = ceil(intProtdistFromCenter);
normDist                  = roundDist / max(roundDist);
[sortNormDist, sortedIdx] = sort(normDist);
sortNormCavInt            = normIntProtInt(sortedIdx);

% Unique values in sorted round distances
uniqueRoundDist = unique(sortNormDist);

% Preallocate arrays
sumIndex        = zeros(1, numel(uniqueRoundDist));
sumIndexIntProt = zeros(1, numel(uniqueRoundDist));

% Total sum of intProt intensities
sumIntProt = sum(sortNormCavInt(:));

% Loop through unique distances
for i = 1:numel(uniqueRoundDist)
    x = uniqueRoundDist(i);
    
    % Find the last index of x in sort_round_dist
    lastIdx     = find(sortNormDist == x, 1, "last");
    sumIndex(i) = lastIdx;
    
    % Sum up intensities up to this index
    sumIndexIntProt(i) = sum(sortNormCavInt(1:lastIdx));
    ratioSumIntProt(i) = sumIndexIntProt(i) / sumIntProt;

end

% Find the index of the closest value
[~, idx25]  = min(abs(ratioSumIntProt - 0.25));
[~, idx50]  = min(abs(ratioSumIntProt - 0.5));
[~, idx75]  = min(abs(ratioSumIntProt - 0.75));
[~, idx100] = min(abs(ratioSumIntProt - 1));

probIntProt = zeros(size(roundScaledIntProtdistFromCent)); % initialize

% Assign probability values based on ranges
probIntProt(roundScaledIntProtdistFromCent > 0     & roundScaledIntProtdistFromCent <= idx25)  = 0.2;
probIntProt(roundScaledIntProtdistFromCent > idx25 & roundScaledIntProtdistFromCent <= idx50)  = 0.4;
probIntProt(roundScaledIntProtdistFromCent > idx50 & roundScaledIntProtdistFromCent <= idx75)  = 0.6;
probIntProt(roundScaledIntProtdistFromCent > idx75 & roundScaledIntProtdistFromCent <= idx100) = 0.8;

% Visualize probIntProt
figure;
imagesc(probIntProt);
%title(['Probability Map of ', proteinName]);
%colorbar;
colormap(customParula);
caxis([0, 1])
axis off;

% Get exactly what MATLAB is displaying:
F      = getframe(gca);      % exact pixels
RGB    = F.cdata;            % exact colors
imwrite(RGB, 'probIntProt.png');  % identical in appearance

%% Section 15: Calculate the density of a given circle for 

% Radial distances (rounded and sorted)
roundDist                  = ceil(intProtdistFromCenter); % Round radial distances
[sortRoundDist, sortedIdx] = sort(roundDist);             % sort roundDist
sortIntProtInt             = normIntProtInt(sortedIdx);   % sort intProt intensity values based on sortRoundDist (keep (x,y) pairs)
uniqueRoundDist            = unique(sortRoundDist);       % Unique radial distances

% Compute annular densities for all unique radial distances
cumulativeDen = zeros(size(uniqueRoundDist)); 

for i = 1:numel(uniqueRoundDist)-1
    rIn  = uniqueRoundDist(i);
    rOut = uniqueRoundDist(i + 1);
    
    area = pi * (rOut^2 - rIn^2); % Corrected annular area formula
    WithinRint = sum(sortIntProtInt(sortRoundDist >= rIn & sortRoundDist <= rOut)); 
    cumulativeDen(i) = WithinRint / area; 
end

% Normalize cumulative density
normCumulativeDen = cumulativeDen / max(cumulativeDen); 

%% **Create Density Heatmap**
denIntProt = zeros(size(roundScaledIntProtdistFromCent)); % Preallocate

% Assign computed densities to corresponding radial regions
for i = 1:numel(uniqueRoundDist)-1
    rIn  = uniqueRoundDist(i);
    rOut = uniqueRoundDist(i + 1);
    
    %if i == 1
       % denIntProt(roundScaledIntProtdistFromCent <= rIn) = normCumulativeDen(i);
    %else
        denIntProt(roundScaledIntProtdistFromCent >= rIn & roundScaledIntProtdistFromCent <= rOut) = normCumulativeDen(i);
    %end
end

%% **Plot Density Heatmap**
figure;
imagesc(denIntProt);
%title(['Density Map of ', proteinName]);
colormap(customHot); % colorbar;
caxis([0, 1]); 
axis off;

% Get exactly what MATLAB is displaying:
F      = getframe(gca);      % exact pixels
RGB    = F.cdata;            % exact colors
imwrite(RGB, 'denIntProt.png');  % identical in appearance

%% Section 18: Save

circleOut.normAvgIntProt = normAvgIntProt;
circleOut.probIntProt    = probIntProt;
circleOut.denIntProt     = denIntProt;
circleOut.outDistVint    = outDistVint;
