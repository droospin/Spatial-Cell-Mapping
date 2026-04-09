
function [data] = four_channel_stacks(source)

%load image file using tiffreadVolume(filename)
imagefile = tiffreadVolume(source);
% imagefile = tiffreadVolume("60x_FM555_CAD488_Phalloidin647_Matrigel_15um_48hr_9_MMStack_Pos0.ome.tif");

%create variable for slices in an image
slices = size(imagefile, 3);

%%

    stackIntProt2 = imagefile( :, :, 3:4:slices);
    intProt2      = max(stackIntProt2, [], 3); 

%%
%create a loop just like above for every 4th slice, starting at 2

    stackRefProt = imagefile( :, :, 4:4:slices);
    refProt      = max(stackRefProt, [], 3);

%%
%create a 3rd loop like 'actin' for every 4th slice, starting at 3

    stackIntProt1 = imagefile( :, :, 2:4:slices);
    intProt1      = max(stackIntProt1, [], 3);

%%
%create a 4th loop like 'cav1' for every 4th slice, starting at 4

    stackDNA = imagefile( :, :, 1:4:slices);
    dna      = max(stackDNA, [], 3);

%% Store variables into a data structure

data.stackIntProt1 = stackIntProt1;
data.stackRefProt  = stackRefProt;
data.stackIntProt2 = stackIntProt2;
data.stackDNA      = stackDNA;
data.intProt1      = intProt1;
data.refProt       = refProt;
data.intProt2      = intProt2;
data.dna           = dna;

end

