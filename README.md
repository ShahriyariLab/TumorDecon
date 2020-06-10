# TumorDecon

Tumors consist of many different cell types, including various immune cells, fibroblasts and epithelial cells. Importantly, each tumor has its own unique variation of immune cells resulting in different responses to the same treatments. Documenting these immune variations will help us to understand the process of tumorigenesis and find ways to arrive at effective treatments.


TumorDecon software includes four deconvolution methods (DeconRNAseq [Gong2013], CIBERSORT [Newman2015], ssGSEA [Şenbabaoğlu2016], Singscore [Foroutan2018]) and several signature matrixes of various cell types, including LM22. The input of this software is the gene expression profile of the tumor, and the output is the relative number of each cell type. Users have an option to choose any of the implemented deconvolution methods and included signature matrixes or import their own signature matrix to get the results.


TumorDecon is available on Github (https://github.com/ShahriyariLab/TumorDecon) and PyPI (https://pypi.org/project/TumorDecon/). To install the software with pip, use both of the following two commands:

```
pip install TumorDecon
pip install git+https://github.com/kristyhoran/singscore
```
