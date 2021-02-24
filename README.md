# TumorDecon


TumorDecon software includes four deconvolution methods (DeconRNAseq [Gong2013], CIBERSORT [Newman2015], ssGSEA [Şenbabaoğlu2016], Singscore [Foroutan2018]) and several signature matrices of various cell types, including LM22. The input of this software is the gene expression profile of the tumor, and the output is the relative number of each cell type. Users have an option to choose any of the implemented deconvolution methods and included signature matrices or import their own signature matrix to get the results. Additionally, TumorDecon can be used to generate customized signature matrices from single-cell RNA-sequence profiles.


TumorDecon is available on Github (https://github.com/ShahriyariLab/TumorDecon) and PyPI (https://pypi.org/project/TumorDecon/). To install the software with pip, use both of the following two commands:

```
pip install TumorDecon
pip install git+https://github.com/kristyhoran/singscore
```

If you use the package or some parts of codes, please cite:
T. Le, R. Aronow, A. Kirshtein, L. Shahriyari, A review of digital cytometry methods: estimating the relative abundance of cell types in a bulk of cells,  Briefing in Bioinformatics, 2020.  
