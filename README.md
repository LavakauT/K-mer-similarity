# K-mer_similarity

K-mer similarity can be applied in important K-mers after MotifDiscovery(https://github.com/ShiuLab/MotifDiscovery) and ML-pipeline (https://github.com/bmmoore43/ML-Pipeline). This repository was developed to replace tamo scripts in the part of "Motif PCC distance using TAMO".

Dap-seq and TFBM motif database of Arabidopsis were regarded as background of similarity comparison. Pearson Correlation Coefficient was the method to determine which kind of known motif is highest similar to each one important K-mers. This repository also included the plot script to generate descriptional heatmap.

> Please prepare the imp files, Dap-seq, and TFBM motif database first before enter any of the scripts. For imp files, there're some sample data you can try in local. Please also notice that database such as Dap-seq included huge amounts of motifs, compare motif between each one K-mer and each motif in database must consumes a lot of resources. We suggested you run this scripts in HPC system.  
