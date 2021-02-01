# scripts:


# renametigs.py
Script to renumber all contigs in an assembly file (identified by ">"-character stating a new line).

# xshift.py
In this script, the spatial euclidean distance between each gene expression profile is calculated. Whereas one genotypes expression profile is artificially temporal shifted. Via iteration, the shifted profile is further shifted since the sum of the distances of all datapoints is on its minimum (assumed artificially perfect correlation). Finally, the applied x-shift is saved. The calculated x-shift does not only enable to track expression shifts on single gene level but can also be compared to phenotypic observed shifts. - Attention: Further improvements to the script are necessary. To this day, each genotype requires well correlated gene expression profiles, just shifted by time.
