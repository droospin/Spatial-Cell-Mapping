
function [data] = three_channel_stacks(source)

%load image file using tiffreadVolume(filename)
imagefile = tiffreadVolume(source);
% imagefile = tiffreadVolume("60x_FM555_CAD488_Phalloidin647_Matrigel_15um_48hr_9_MMStack_Pos0.ome.tif");

%create variable for slices in an image
slices = size(imagefile,3);

%%

    stackIntProt = imagefile(:,:,2:3:slices);
    intProt      = max(stackIntProt, [], 3);

%%
    stackRefProt = imagefile(:,:,3:3:slices);
    refProt      = max(stackRefProt, [], 3);

%%
    stackDNA = imagefile(:,:,1:3:slices);
    dna      = max(stackDNA, [], 3);

%%
data.stackDNA     = stackDNA;
data.stackRefProt = stackRefProt;
data.stackIntProt = stackIntProt;
data.dna          = dna;
data.intProt      = intProt;
data.refProt      = refProt;

end