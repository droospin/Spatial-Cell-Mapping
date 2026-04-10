# Spatial-Cell-Mapping
MATLAB & R-based computational mapping tools used to characterize and map cells confined to micropattern geometries that are used to force cells into different physiological behaviors like migration and stabilization.

The power of micropatterns comes from the idea that all cells being imaged will have the same shape, so you can easily align cells and look at trends in protein localizations and structural organizations during specific behaviors, something that is much more difficult to do with most standardly used cell culture techniques.

For more information, check out the following:

Grespin, D. B., Farrington, J. S., Niven, T. G., Russell, L. J., Loerke, D., David, A. J., Grespin, M. S., Culkin, C. M., Bartoletti, A. P., Meadows, S. M., & Kushner, E. J. (2026). PRINCIPLES GOVERNING ENDOTHELIAL CAVEOLAE ORGANIZATION DURING ANGIOGENESIS. bioRxiv : the preprint server for biology, 2026.03.27.714916. https://doi.org/10.64898/2026.03.27.714916

Grespin, D. B., Niven, T. G., Babson, R. O., & Kushner, E. J. (2023). Lipidure-based micropattern fabrication for stereotyping cell geometry. Scientific reports, 13(1), 20451. https://doi.org/10.1038/s41598-023-47516-8

The stack readers are meant for reading z-stack tiff files in MATLAB. If you need to read ND2 image files, I would considering using this very well organized Bioformats Image Toolbox assembled by Jian Wei Tay: https://matlab.mathworks.com/open/fileexchange/v1?id=129249 

The functions in the main branch are important for most of the functions in other branches to run. What I do, due to their ubiquitous nature and to not add too much code to one function, is save them in my path that I am operating in for the function. Then before I run the main function, I ensure that those files are on my working directory or path and the main function.

Please reach out with questions about example datasets, using MATLAB or R, and/or about the code. Cheers!
