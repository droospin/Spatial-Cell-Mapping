function [result] = smooth_and_plot(input, alpha, titleStg, xAxisStg, yAxisStg, currXlim, currYlim, Xtick, Ytick)
% smooth_and_plot: Plots an interpolated version of x and y data, with a standard deviation
% error envelope (can be adjusted to plot 'standard error of the mean'
% error envelope instead, if desired)
%
% INPUT:    input = n x 2 matrix with x (column 1) and y (column 2) data to
%                   be plotted.
% 
% OUPUT:   result = n x 2 matrix containing:
%                   first column:   binned x-coordinates
%                   second column:  interpolated y-coordinates
%                   third column:   standard deviation
%                   fourth column:  n (sum of normalized weights)
%                   fifth column:   standard error of mean (=std/sqrt(n))
%          Output plot also generated, with both raw and interpolated versions.
%
% L. Russell w/ minor modification from D. Grespin (updated March 2025)
od = cd; % original directory

%% CONSTANTS
% figure formatting strings (change these to match what you're plotting):

% % rgb triple for plotting color (uncomment desired color):
% rgb_vec = [1 0 0]; % red
% rgb_vec = [0 1 0]; % green
rgbVec = [0 0 1]; % blue
% rgb_vec = [1 0 1]; % magenta
% rgb_vec = [0 1 1]; % cyan
% rgb_vec = [1 1 0]; % yellow
% rgb_vec = [0.5 0.5 0.5]; % grey

%% MAIN
num_bins = 50; % # of equally sized bins for interpolation, change as needed

% find range of x-data, to select appropriately sized bins
range_x = range(input(:,1)); % x-range
x_step = range_x / num_bins; % divide range up into num_bins equal steps
x_bins = min(input(:,1)):x_step:max(input(:,1)); % bin vector, for interpolation

sigma = 2 * x_step; % smoothing factor, change as needed (higher value = more smoothing)

% perform interpolation:
[result] = interpolatePointVals2Vec_gauss_wstd(input, x_bins, sigma, 0);

%% PLOTTING
figure; % initialize figure object

% scatter plot of raw data
scatter(input(:,1), input(:,2), 2, 'k', 'filled', 'MarkerFaceAlpha', alpha);
hold on;
% figure aesthetics/ labels
%set(gca, 'FontName', 'Arial', 'FontSize', 17);
% title(titleStg, 'Fontweight', 'normal', 'FontName', 'Arial', 'FontSize', 20);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
ax.LineWidth = 1; % make the grid thicker (default is 0.5)
grid off;
axis tight; box off;
% xlabel(xAxisStg, 'FontName', 'Arial', 'FontSize', 18);
% ylabel(yAxisStg, 'FontName', 'Arial', 'FontSize', 18);
xlim(currXlim); % xticks(Xtick);
ylim(currYlim); % yticks(Ytick);
fig = gcf;
fig.Units = 'pixels';
fig.Position = [0, 0, 264, 200];  % [left, bottom, width, height]

figure; % initialize figure object

% plot interpolated line:
plot(result(:,1), result(:,2), 'Color',rgbVec); hold on; grid on;

% plot error envelope (2 options):
% Standard deviation (this is what I generally use, unless data is particularly noisy)
errorenvelope(result(:,1)',result(:,2)', result(:,3)' ./ 2, rgbVec, 0.3);
% Standard error of the mean (uncomment below if desired):
% errorenvelope(result(:,1)',result(:,2)', result(:,5)', rgb_vec, 0.3);

% figure aesthetics/ labels
%set(gca, 'FontName', 'Arial', 'FontSize', 17);
% title({titleStg,'(interpolated)'}, 'FontWeight', 'normal', 'FontName', 'Arial', 'FontSize', 10);
ax.LineWidth = 1; % make the grid thicker (default is 0.5)
grid on;
axis tight; box off;
% xlabel(xAxisStg, 'FontName', 'Arial', 'FontSize', 9); 
% ylabel(yAxisStg, 'FontName', 'Arial', 'FontSize', 9);
xlim(currXlim); xticks(Xtick);
ylim(currYlim); yticks(Ytick);
fig = gcf;
fig.Units = 'pixels';
fig.Position = [0, 0, 264, 200];  % [left, bottom, width, height]

cd(od); % return to original directory

end % end smooth_and_plot function


%% HELPER FUNCTIONS:

function [ result ] = interpolatePointVals2Vec_gauss_wstd( valuematrix, bins, sigma, displayVariable)
% construct a value vector using values at unevenly spaced locations; the 
% parameter vector is interpolated using gaussian weighting of distances; 
% in addition it is possible to weight contributions using an external
% weight vector
%
% INPUT:    valuematrix     n x 2 matrix containing:
%                           first column:   features' x-coordinates
%                           second column:  features' y-coordinates
%                           third column: (optional) weight
%
%           bins            vector of the desired value bins of the
%                           x-vector to which the value will be
%                           interpolated
%
%           sigma           sigma of the Gaussian weighting function (pix)
%                           smaller sigmas will of course allow you to 
%                           resolve smaller spatial variations, but it
%                           makes no sense to choose a sigma that is much
%                           below the typical distance between individual
%                           features
%
% OUPUT:   result           n x 2 matrix containing:
%                           first column:   bin x-coordinates
%                           second column:  interpolated y-coordinates
%                           third column:   standard deviation
%                           fourth column:  n (sum of normalized weights)
%                           fifth column:   standard error of mean (=std/sqrt(n))
%
% written 02/06/2014 Dinah Loerke
% updated 08/21/2020 DL for n and sem

display = 0;
if nargin>3
    if displayVariable==1
       display = 1;
    end
end

%   extract x,y coordinate and parameter vectors from input
xvec = valuematrix(:,1);
yvec = valuematrix(:,2);

weightvec_sig = 1+0*xvec;
if size(valuematrix,2)>2
   weightvec_sig = valuematrix(:,3);
end

% initialize results
result = nan*zeros(length(bins),5);

% loop over bins
for i=1:length(bins)
    
    % display progress
    % display(['processing mask section ',num2str(i),' of ', num2str(length(stepvector))]);
      
    % distance vector (from current bin) for each data point
    distvec = abs(xvec-bins(i));
    
    % weighting vector for each data point (Gaussian function with distance
    % from bin)
    weightvec_dist = exp(-(distvec.^2)/(2*sigma^2));
    
    nanpos = isnan(yvec);
    weightvec_dist(nanpos) = nan;
    
    % multiply this weight vector with the 'external' weights, which may
    % represent e.g. statistical significance of the point
    weightvec = weightvec_dist .* weightvec_sig;
    
    % weighted average of the data parameter (yvec) for this bin
    valuvec = nansum((yvec.*weightvec),1)./nansum(weightvec,1);
    
    % to estimate standard deviation, first determine square difference of
    % each point from average 
    varvec = (yvec - valuvec).^2;
    % to estimate variance (and then std), again weight these differences
    % with the weightvector
    variance = nansum((varvec.*weightvec),1)./nansum(weightvec,1);
    
    % n = sum of normalized weights
    num = sum(weightvec,1, 'omitnan');
     
    % write results (averaged parameter values for the current pixels) into
    % the results at the specified positions 
    result(i,1:5) = [bins(i) valuvec sqrt(variance) num sqrt(variance)/sqrt(num)];
    
end % of for i-loop

% display results
if display == 1
    figure;
    hold off; plot(xvec,yvec,'bo');
    hold on; plot(result(:,1),result(:,2),'r.-');
    hold on; errorbar(result(:,1),result(:,2),result(:,3),'g-');
    hold on; errorbar(result(:,1),result(:,2),result(:,5),'m-');
end

end % end interpolation function

%%

function errorenvelope(x,y,err,faceColor,alphaVal)
% errorenvelope generates a symmetric error envelope, y-err to y+err, over
% the user-input x-vector.
%
% INPUTS:               x: vector of x-values
%                       y: vector of y-values which is the mean of the
%                          envelope
%                     err: vector of error values which will be half the
%                          width of the envelope ( y +/- err).
%               faceColor: optional input of the envelope color, dfault is
%                          blue
%                alphaVal: transparency value where 1 is opaque and zero is
%                          transparent, the default is 0.2;
%
%

% if no facecolor is given as input, set to blue
if nargin < 4
    faceColor = [0 0 1];
end
if nargin < 5
    alphaVal = 0.2;
end

% we use an error envelope that is symmetric about the mean
ymin = y - err;
ymax = y + err;

% if NaN's exist in the y-values, take just the finite values
finitePts  = isfinite(ymin);
xfinite    = x(finitePts);
yminfinite = ymin(finitePts);
ymaxfinite = ymax(finitePts);

x    = xfinite;
ymin = yminfinite;
ymax = ymaxfinite;

X = [x,fliplr(x)];                % create continuous x value array for plotting
Y = [ymin,fliplr(ymax)];          % create y values for out and then back
h = fill(X,Y,'b');   
set(h, 'FaceColor', faceColor, 'EdgeColor', 'none')

% make the filled area transparent, this can be adjusted by changing the
% value: 1 is opaque, and 0 is completely transparent
alpha(alphaVal)

end % end error envelope function