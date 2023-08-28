# NetID
# Introduction
Liquid chromatography–high-resolution mass spectrometry (LC-MS)-based metabolomics aims to identify and quantify all metabolites, but most LC-MS peaks remain unidentified. Here we present a global network optimization approach, NetID, to annotate untargeted LC-MS metabolomics data. The approach aims to generate, for all experimentally observed ion peaks, annotations that match the measured masses, retention times and (when available) tandem mass spectrometry fragmentation patterns. Peaks are connected based on mass differences reflecting adduction, fragmentation, isotopes, or feasible biochemical transformations. Global optimization generates a single network linking most observed ion peaks, enhances peak assignment accuracy, and produces chemically informative peak–peak relationships, including for peaks lacking tandem mass spectrometry spectra. Thus, NetID applies existing metabolomic knowledge and global optimization to substantially improve annotation coverage and accuracy in untargeted metabolomics datasets, facilitating metabolite discovery.

NetID requires: (1) data file (in .mzXML format), (2) a peak table (in .csv format), (3) a reference compound library (in .rds format), (4) a transformation table (in .csv format), for which we assembled a list of 25 biochemical atom differences and 59 abiotic atom differences. NetID optionally use (5) a list of known metabolites' retention time, for which we provide our in-house retention time list for demonstration and (6) .mgf file containing MS2 information and (7) MS2 reference library (in .rds format). More details in section 3.1. 

In the 2023 August version, we implemented the following updates:  
(1) Include both orbitrap and TOF demo data;  
(2) Handle data-dependent MS2 data;  
(3) Incorporate FastNetID codes proposed in OmicsNet<sup>1</sup>, which uses C language to accelerate calculations in NetID;  
(4) Provide a detailed peak picking workflow using MZmine3<sup>2</sup> and EVA<sup>3</sup>;  
(5) Merge two compound libraries HMDB<sup>4</sup> and PubChemLite<sup>5-6</sup>;  
(6) Restructure the workflow to remove unstable functions and simplify the parameter settings.

**Citation**: Chen, L., Lu, W., Wang, L. et al. Metabolite discovery through global annotation of untargeted metabolomics data. Nat Methods 18, 1377–1385 (2021). https://doi.org/10.1038/s41592-021-01303-3  
**Git-hub**: https://github.com/LiChenPU/NetID

**Reference**: 

1. Zhou, G., Pang, Z., Lu, Y., Ewald, J. & Xia, J. OmicsNet 2.0: A web-based platform for multi-omics integration and network visual analytics. Nucleic Acids Research gkac376 (2022).
2. Schmid, R. et al. Integrative analysis of multimodal mass spectrometry data in MZmine 3. Nat Biotechnol 41, 447–449 (2023).
3. Guo, J. et al. EVA: Evaluation of Metabolic Feature Fidelity Using a Deep Learning Model Trained With Over 25000 Extracted Ion Chromatograms. Anal. Chem. 93, 12181–12186 (2021).
4. Wishart, D. S. et al. HMDB 4.0: The human metabolome database for 2018. Nucleic Acids Research 46, D608–D617 (2018).
5. Kim, S. et al. PubChem 2019 update: Improved access to chemical data. Nucleic Acids Research 47, D1102–D1109 (2019).
6. Bolton, E. & Schymanski, E. PubChemLite tier0 and tier1. Zenodo https://doi.org/10.5281/zenodo.3611238 (2020).
