
%% For stable vessels

for i = 1:length(mouseStableLocker)

    refProt     = mouseStableLocker(i).refProt;
    refProtEC   = mouseStableLocker(i).refProtEC;
    intProt     = mouseStableLocker(i).intProt;
    intProtEC   = mouseStableLocker(i).intProtEC;

    normRefProt = double(refProt) / double(max(refProt(:)));
    normIntProt = double(intProt) / double(max(intProt(:)));

    normRefProtEC                 = normRefProt;
    normRefProtEC(refProtEC == 0) = 0;

    normIntProtEC                 = normIntProt;
    normIntProtEC(refProtEC == 0) = 0;

    mouseStableInt(i).normRefProtEC = normRefProtEC;
    mouseStableInt(i).normIntProtEC = normIntProtEC;

end

refProtECvec   = [];
intProtECvec   = [];

for i = 1:numel(mouseStableInt)
           
    refProtECvec = [refProtECvec; mouseStableInt(i).normRefProtEC(:)];
    intProtECvec = [intProtECvec; mouseStableInt(i).normIntProtEC(:)];
    idx          = refProtECvec > 0 & refProtECvec <= 1;

end

mouseStableInt(1).refProtVintProt = [refProtECvec(idx), intProtECvec(idx)];

titleStg       = 'Normalized Vecad vs. Cav1 Intensity';
xAxisStg       = 'VE-cad';
yAxisStg       = 'Cav1';
lineXlim       = [0 1]; Xtick = 0:0.2:1;
lineYlim       = [0 1]; Ytick = 0:0.2:1;

[mouseStableInt(1).outRefProtVintProt] = smooth_and_plot(mouseStableInt(1).refProtVintProt, 0.02, titleStg, xAxisStg, yAxisStg, lineXlim, lineYlim, Xtick, Ytick);

figure; histogram(mouseStableInt(1).refProtVintProt(:, 1), 10)
figure; histogram(mouseStableInt(1).refProtVintProt(:, 2), 10)
%% For tip cells

for i = 1:length(mouseTipLocker)

    refProt     = mouseTipLocker(i).refProt;
    refProtEC   = mouseTipLocker(i).refProtEC;
    refProtDist = mouseTipLocker(i).refProtDist;
    intProt     = mouseTipLocker(i).intProt;
    intProtEC   = mouseTipLocker(i).intProtEC;
    intProtDist = mouseTipLocker(i).intProtDist;

    normRefProt     = double(refProt) / double(max(refProt(:)));
    normRefProtDist = refProtDist     / max(refProtDist(:));
    normIntProt     = double(intProt) / double(max(intProt(:)));
    normIntProtDist = intProtDist     / max(intProtDist(:));

    normRefProtEC                 = normRefProt;
    normRefProtEC(refProtEC == 0) = 0;

    normIntProtEC                 = normIntProt;
    normIntProtEC(refProtEC == 0) = 0;

    mouseTipOut(i).normRefProtEC   = normRefProtEC;
    mouseTipOut(i).normRefProtDist = normRefProtDist;
    mouseTipOut(i).normIntProtEC   = normIntProtEC;
    mouseTipOut(i).normIntProtDist = normIntProtDist;

end

refProtECvec   = [];
refProtDistVec = [];
intProtECvec   = [];
intProtDistVec = [];

for i = 1:numel(mouseTipOut)
           
    refProtECvec   = [refProtECvec;   mouseTipOut(i).normRefProtEC(:)];
    refProtDistVec = [refProtDistVec; mouseTipOut(i).normRefProtDist(:)]; 
    intProtECvec   = [intProtECvec;   mouseTipOut(i).normIntProtEC(:)];
    intProtDistVec = [intProtDistVec; mouseTipOut(i).normIntProtDist(:)];
    idx            = refProtECvec > 0 & refProtECvec <= 1;

end

mouseTipOut(1).refProtVintProt = [refProtECvec(idx), intProtECvec(idx)];

titleStg       = 'Normalized Vecad vs. Cav1 Intensity';
xAxisStg       = 'VE-cad';
yAxisStg       = 'Cav1';
lineXlim       = [0 1]; Xtick = 0:0.2:1;
lineYlim       = [0 1]; Ytick = 0:0.2:1;

[mouseTipOut(1).outRefProtVintProt] = smooth_and_plot(mouseTipOut(1).refProtVintProt, 0.02, titleStg, xAxisStg, yAxisStg, lineXlim, lineYlim, Xtick, Ytick);

mouseTipOut(1).intProtDistVInt = [intProtDistVec(idx), intProtECvec(idx)];

figure; histogram(mouseTipOut(1).refProtVintProt(:, 1), 10)
figure; histogram(mouseTipOut(1).refProtVintProt(:, 2), 10)

titleStg       = 'Cav1 Normalized Distance vs. Intensity';
xAxisStg       = 'Distance from Migratory Front';
yAxisStg       = 'Cav1 Intensity';
lineXlim       = [0 1]; Xtick = 0:0.2:1;
lineYlim       = [0 1]; Ytick = 0:0.2:1;

[mouseTipOut(1).outIntProtDistVInt] = smooth_and_plot(mouseTipOut(1).intProtDistVInt, 0.02, titleStg, xAxisStg, yAxisStg, lineXlim, lineYlim, Xtick, Ytick);

figure; histogram(mouseTipOut(1).intProtDistVInt(:, 1), 10)
figure; histogram(mouseTipOut(1).intProtDistVInt(:, 2), 10)