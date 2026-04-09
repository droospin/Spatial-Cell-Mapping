function [vesselMask, outDistVInt] = read_retinal_nd2(folder)

    vesselIm    = process_nd2_folder(folder);
    vesselMask  = cellMaskingGUI(vesselIm);  % Only get fig first

end

%%
function vesselData = process_nd2_folder(folder)

    imFiles = fullfile(folder, '**', '*.nd2');
    fileList = dir(imFiles);

    for i = 1:length(fileList)
        source = fullfile(fileList(i).folder, fileList(i).name);
        imObj  = BioformatsImage(source);

        xSize  = imObj.height;
        ySize  = imObj.width;

        dna     = getPlane(imObj, 1, 'Cy5', 1);
        refProt = getPlane(imObj, 1, 'FITC', 1);
        intProt = getPlane(imObj, 1, 'TRITC', 1);

        % Display
        %figure; imshowpair(intProt, refProt); pause(0.5);
        %figure; imshowpair(intProt, refProt, 'montage'); pause(0.5);
        
        % close displays
        close all;

        % Store
        vesselData(i).filename     = fileList(i).name;
        vesselData(i).intProt      = intProt;
        vesselData(i).refProt      = refProt;
        vesselData(i).dna          = dna;
        vesselData(i).xSize        = xSize;
        vesselData(i).ySize        = ySize;
    end

end

%%

function vesselMask = cellMaskingGUI(inputData)

    vesselMask = inputData;  % Copy all fields

    for n = 1:numel(vesselMask)
        vesselMask(n).mask_indices = {};
        vesselMask(n).intProtSeg   = [];
    end

    currentIdx = 1;
    fig = figure('Name', 'Cell Masking GUI', 'Position', [100, 100, 1400, 1200], ...
                 'MenuBar', 'none', 'ToolBar', 'none');
    ax = axes('Parent', fig, 'Units', 'pixels', 'Position', [100, 0, 1024, 1024]);

    btnNext = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'Next', ...
                    'Position', [1200, 750, 80, 30], 'Callback', @(src, event) navigate(1));
    btnPrev = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'Previous', ...
                    'Position', [1290, 750, 80, 30], 'Callback', @(src, event) navigate(-1));

    btnTrace = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'Trace Mask', ...
                     'Position', [1200, 700, 170, 30], 'Callback', @(src, event) trace_and_mask());

    btnDone = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'Done', ...
                    'Position', [1200, 650, 170, 30], 'Callback', @(src, event) closeGUI());

    btnClear = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'Clear Masks', ...
                     'Position', [1200, 600, 170, 30], 'Callback', @(src, event) clear_masks());

    mask = zeros(size(inputData(1).refProt));
    updateDisplay();
    
    uiwait(fig);  % This will pause execution after the display is ready


    % --- NESTED FUNCTIONS ---
    function updateDisplay()
        cla(ax);
        dna     = inputData(currentIdx).dna;
        refProt = inputData(currentIdx).refProt;
        intProt = inputData(currentIdx).intProt;
    
        % Normalize and boost refProt
        refProtNorm = mat2gray(refProt);
        intProtNorm = mat2gray(intProt);
        dnaNorm     = mat2gray(dna);
        
        boost = 1.5;  % Try between 1.1 and 2.0
        dnaAlpha = 0.5;  % Reduce contribution of intProt
        
        rgb_image = cat(3, ...
                        min(refProtNorm * boost, 1), ... % Red channel (refProt)
                        zeros(size(intProtNorm)), ... % Green channel (intProt dimmed)
                        min(dnaNorm * dnaAlpha, 1));     % Blue channel (not used)
    
        imshow(rgb_image, [], 'Parent', ax);
        
        % Overlay current mask
        if ~isempty(vesselMask(currentIdx).mask_indices)
            hold(ax, 'on');
            for i = 1:numel(vesselMask(currentIdx).mask_indices)
                boundary = vesselMask(currentIdx).mask_indices{i};
                plot(ax, boundary(:,2), boundary(:,1), 'y', 'LineWidth', 1.5);
            end
            hold(ax, 'off');
        end
    end


    function navigate(step)
        currentIdx = mod(currentIdx - 1 + step, length(inputData)) + 1;
        updateDisplay();
    end


    function trace_and_mask()
        h = drawfreehand(ax, 'Color', 'y');
        new_mask = createMask(h);
        delete(h);
    
        save_curr_mask(new_mask);  % Pass it in here
        updateDisplay();
    end
    
    function save_curr_mask(new_mask)
        % Add to current list
        B = bwboundaries(new_mask);
        if isempty(B)
            return;
        end
    
        % Combine old and new
        vesselMask(currentIdx).mask_indices = ...
            [vesselMask(currentIdx).mask_indices; B];
    
        % Generate full mask
        cumulative_mask = false(size(new_mask));
        for k = 1:numel(vesselMask(currentIdx).mask_indices)
            boundary = vesselMask(currentIdx).mask_indices{k};
            cumulative_mask = cumulative_mask | ...
                poly2mask(boundary(:,2), boundary(:,1), size(new_mask,1), size(new_mask,2));
        end
    
        % Save new segmented intensity
        intProt = inputData(currentIdx).intProt;
        vesselMask(currentIdx).intProtSeg = intProt .* cast(cumulative_mask, class(intProt));
        vesselMask(currentIdx).binaryMask = cumulative_mask;

    end

    function clear_masks()
        mask = zeros(size(inputData(currentIdx).refProt));
        vesselMask(currentIdx).mask_indices = {};
        vesselMask(currentIdx).intProtSeg   = {};
        updateDisplay();
    end

    % In closeGUI:
    function closeGUI()
        assignin('base', 'stableVessel', vesselMask);
        uiresume(fig);  % Resume main script
        delete(fig);    % Close GUI
    end

end