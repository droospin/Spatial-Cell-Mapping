# Confluence-Track-and-Lines
Please read the main readme file before using this or any other pipeline.

This pipeline is built for confluency experiments, particularly micropatterns. 

This is pipeline is built to analyze and map monolayer cell z-stacks on tracks and line micropatterns with four channels: a reference channel (most often F-actin, to mark the cell outlines/junctions), a junctional a protein of interest, a protein of interest, and DNA (Hoechst dye).

To ensure you are analyzing the desired channel, use the proper channel reader in the main branch and align the numbers with the channels you see.

There are two pipelines in this branch: 1) line_confluency_mapping can be used for standard pipeline tracing individual cells in both track patterns and line patterns 2) line_confluency_junction_mapping will allow you to trace the junctions across the entire image (more of a global than a per cell quanitification).
