# Mouse-Retinal-Vasculature
Please read the main readme file before using this or any other pipeline.

There are two main pipelines for this branch. Additionally, the functions in this branch (and this branch only) should be self-sufficient to run the code.

The first main pipeline is mouseModelAnalysis. This function will guide you to trace a line across the vascular/angiogenic/migratory front. Then it will calculate the intensities of the protein of interest as a function of normalized distance away from the front. 

The second main pipline mouse_per_cell. As it sounds, this function will include tracing cells close to the angiogenic front and mapping the intensities against corresponding distances from the angiogenic front. It is designed to help you see whether a protein/image signal is more rear or front localized in individual cells.
