# Line-Patterns-Single
Please read the main readme file before using this or any other pipeline.

This is pipeline is built to analyze and map individual cell z-stacks on line micropatterns with three channels: a reference channel (most often F-actin or a membrane stain, to mark the cell shape), a protein of interest, and DNA (Hoechst dye).

To ensure you are analyzing the desired channel, use the proper channel reader in the main branch and align the numbers with the channels you see.

There are two pipelines in this branch: 1) line_alignment is meant to be run on an image folder, with its output being the input for line_analysis 2) line_3d as a stand alone with expected input of an image folder.
