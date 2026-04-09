function mouseTipLast = mouse_per_cell(mouseTipLocker, mouseTipOut)
%% Prepare combined images per sample
v1 = mouseTipLocker;
v2 = mouseTipOut;
mouseTipLast = struct([]);

for i = 1:length(v1)
    dna = double(v1(i).dna);
    normDNA = dna / max(dna(:));

    mouseTipLast(i).filename      = v1(i).filename;
    mouseTipLast(i).normDNA       = normDNA;
    mouseTipLast(i).normRefProtEC = v2(i).normRefProtEC;
    mouseTipLast(i).normIntProtEC = v2(i).normIntProtEC;

    % Combine into RGB for visualization
    mouseTipLast(i).comboIm = cat(3, ...
        mouseTipLast(i).normRefProtEC, ...
        mouseTipLast(i).normIntProtEC, ...
        mouseTipLast(i).normDNA);
end

%% Optional: Quick overview of first 10 images
figure;
tiledlayout(2,5,'TileSpacing','tight','Padding','tight');
for i = 1:min(10,length(mouseTipLast))
    nexttile
    imshow(mouseTipLast(i).comboIm,[]);
    title(sprintf('Image %d', i))
end
pause(2); close;

mouseTipLast(1).intProtTotNormDistVInt = [];
mouseTipLast(1).refProtTotNormDistVInt = [];

%%
for i = 1:length(mouseTipLast)
    combo = mouseTipLast(i).comboIm;
    intensityImRefProt = combo(:,:,1);  % adjust channel if needed
    intensityImIntProt = combo(:,:,2);

    intProtDistVec = [];
    intProtIntVec  = [];

    refProtDistVec = [];
    refProtIntVec  = [];

    normDist   = [];
    intProtInt = [];
    refProtInt = [];

    hFig = figure('Name','Per-Cell Front Analysis','NumberTitle','off');
    hAx  = axes('Parent',hFig);
    imshow(combo,[],'Parent',hAx);
    hold(hAx,'on');

    doneImage = false;  % <--- new: control flag

    % Create buttons
    uicontrol('Style','pushbutton','String','Draw Front + Cell',...
        'Position',[20 20 140 40],...
        'Callback',@drawFrontThenCell);

    uicontrol('Style','pushbutton','String','Done Image',...
        'Position',[180 20 100 40],...
        'Callback',@(src,evt) setappdata(hFig,'doneImage',true));

    % --- NEW: keep figure open until user clicks Done Image ---
    setappdata(hFig,'doneImage',false);
    while ~getappdata(hFig,'doneImage')
        pause(0.1);  % allow GUI callbacks to run
    end
    close(hFig);

    mouseTipLast(i).intProtNormDistVInt = [intProtDistVec intProtIntVec];
    mouseTipLast(i).refProtNormDistVInt = [refProtDistVec refProtIntVec];

    mouseTipLast(1).intProtTotNormDistVInt = [mouseTipLast(1).intProtTotNormDistVInt; mouseTipLast(i).intProtNormDistVInt];
    mouseTipLast(1).refProtTotNormDistVInt = [mouseTipLast(1).refProtTotNormDistVInt; mouseTipLast(i).refProtNormDistVInt];
    
    mouseTipLast(i).normDist   = normDist;
    mouseTipLast(i).intProtInt = intProtInt;
    mouseTipLast(i).refProtInt = refProtInt;

end

mouseTipLast(1).intProtInterp = smooth_and_plot(mouseTipLast(1).intProtTotNormDistVInt, 0.25, '', '', '', [0 1], [0 1], 0:0.25:1, 0:0.25:1);
mouseTipLast(1).refProtInterp = smooth_and_plot(mouseTipLast(1).refProtTotNormDistVInt, 0.25, '', '', '', [0 1], [0 1], 0:0.25:1, 0:0.25:1);

%% -------- Nested function: Draw Front + Cell --------
function drawFrontThenCell(~,~)

    % Draw Front
    title(hAx,'Draw migratory front. Double-click to finish.');
    hLine         = drawpolyline('Parent',hAx);
    frontVertices = hLine.Position;

    canvas = zeros(size(combo, 1),size(combo, 2),'uint8');
    for n = 1:size(frontVertices,1)-1
        canvas = insertShape(canvas,'Line',...
                 [frontVertices(n,:) frontVertices(n+1,:)],...
                 'Color','white','LineWidth',2);
    end

    frontMask = rgb2gray(canvas) > 0;
    distMap   = bwdist(frontMask); distMap(frontMask) = 0;
    plot(hAx,frontVertices(:,1),frontVertices(:,2),'w','LineWidth',2);

    % Draw Cell
    title(hAx,'Draw ONE cell polygon.');
    hCell    = drawpolygon('Parent',hAx);
    cellMask = createMask(hCell);

    cellDist = distMap(cellMask);
    if isempty(cellDist), return; end

    % --- Normalize distances per cell ---
    dmin = min(cellDist); 
    dmax = max(cellDist);
    if dmax==dmin
       normDist = zeros(size(cellDist));  % all same distance
    else
       normDist = (cellDist - dmin) ./ (dmax - dmin);  % 0 to 1 per cell
    end
        
    % --- Normalize IntProt intensity per cell ---
    intProtInt = intensityImIntProt(cellMask);
    maxInt     = max(intProtInt(:));
    if maxInt == 0
       normIntIntProt = zeros(size(intProtInt));
    else
       normIntIntProt = intProtInt / maxInt;  % per-cell max = 1
    end

    % --- Normalize RefProt intensity per cell ---
    refProtInt = intensityImRefProt(cellMask);
    maxInt     = max(refProtInt(:));
    if maxInt == 0
       normIntRefProt = zeros(size(refProtInt));
    else
       normIntRefProt = refProtInt / maxInt;  % per-cell max = 1
    end
    
    % Only keep pixels > 0 (optional, to match previous logic)
    idx = normIntIntProt > 0;
        
    % Append normalized values to vectors
    intProtDistVec = [intProtDistVec; normDist(idx)];
    intProtIntVec  = [intProtIntVec;  normIntIntProt(idx)];

    refProtDistVec = [refProtDistVec; normDist(idx)];
    refProtIntVec  = [refProtIntVec;  normIntRefProt(idx)];

    % Visual feedback
    plot(hAx,hCell.Position(:,1),hCell.Position(:,2),'y','LineWidth',1.5);
    title(hAx,sprintf('Cells saved: %d', size(intProtDistVec,1)));
    end
end