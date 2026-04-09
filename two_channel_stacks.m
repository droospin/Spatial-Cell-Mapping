
function [data] = two_channel_stacks(source)

%load image file using tiffreadVolume(filename)
imagefile = tiffreadVolume(source);

%create variable for slices in an image
slices = size(imagefile,3);
stackRefProt = zeros(512,512,ceil(slices/2));
stackIntProt = zeros(512,512,ceil(slices/2));
%%
%create a loop that collects every 2nd slice and places them into 'dna',
%a stack that is every 2nd slice at starting at 1. 

%ceil rounds up to nearest whole number

for iref = 1:2:slices
    stackRefProt = imagefile(:,:,1:2:iref);
    refProt = max(stackRefProt, [], 3);                % maximum z-projection of stack_dna
end
%%
%create a loop just like above for every 2nd slice, starting at 2

for iint = 2:2:slices
    stackIntProt = imagefile(:,:,2:2:iint);
    intProt = max(stackIntProt, [], 3);            % maximum z-projection of stack_actin
end

%%

data.stackRefProt = stackRefProt;
data.stackIntProt = stackIntProt;
data.refProt = refProt;
data.intProt = intProt;

end

