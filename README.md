# PhD thesis: statistical methods for the integrative analysis of single-cell multi-omics data

- Author: Ricard Argelaguet  
- Supervisors: John Marioni and Oliver Stegle  

<center>
<img src="Figs/both_logos.png" width="400"/>
</center>

## Download

- [University of Cambridge Apollo's repository](https://www.repository.cam.ac.uk/handle/1810/315822)  
- [Dropbox link](https://www.dropbox.com/s/4j7gasasvck206j/thesis.pdf?dl=0)


## Introduction

Single-cell profiling techniques have provided an unprecedented opportunity to study cellular heterogeneity at the molecular level. This represents a remarkable advance over traditional bulk sequencing methods, particularly to study lineage diversification and cell fate commitment events in heterogeneous biological processes. While the large majority of single-cell studies are focused on quantifying RNA expression, transcriptomic readouts provide only a single dimension of cellular heterogeneity. Recently, technological advances have enabled multiple biological layers to be probed in parallel one cell at a time, unveiling a powerful approach for investigating multiple dimensions of cellular heterogeneity. However, the increasing availability of multi-modal data sets needs to be accompanied by the development of suitable integrative strategies to fully exploit the data generated. In this thesis I worked in collaboration with different research groups to introduce innovative experimental and computational strategies for the integrative study of multi-omics at single-cell resolution.

## Contributions

**The first contribution is the development of scNMT-seq, a protocol for the simultaneous profiling of RNA expression, DNA methylation and chromatin accessibility in single cells**. We demonstrate how this assay provides a powerful approach for investigating regulatory relationships between the epigenome and the transcriptome within individual cells. [Paper](https://www.nature.com/articles/s41467-018-03149-4) and [github repository (outdated)](https://github.com/PMBio/scNMT-seq)

**The second contribution is Multi-Omics Factor Analysis (MOFA), a statistical framework for the unsupervised integration of multi-omics data sets**. MOFA is a Bayesian latent variable model that can be viewed as a statistically rigorous generalization of Principal Component Analysis to multi-omics data. The method provides a principled approach to retrieve, in an unsupervised manner, the underlying sources of sample heterogeneity while at the same time disentangling which axes of heterogeneity are shared across multiple modalities and which are specific to individual data modalities. [Paper](https://www.embopress.org/doi/10.15252/msb.20178124) and [webpage](https://biofam.github.io/MOFA2/)

**The third contribution is the generation of a comprehensive molecular roadmap of mouse gastrulation at single-cell resolution using scNMT-seq**. We employed scNMT-seq to simultaneously profile RNA expression, DNA methylation and chromatin accessibility for hundreds of cells, spanning multiple time points from the exit from pluripotency to primary germ layer specification. Using MOFA, and other tools, we performed an integrative analysis of the multi-modal measurements, revealing novel insights into the role of the epigenome in regulating this key developmental process. [Paper](https://www.nature.com/articles/s41586-019-1825-8) and [github repository](https://github.com/rargelaguet/scnmt_gastrulation)

**The fourth contribution is MOFA+, an extended formulation of the MOFA model tailored to the analysis of large-scale single-cell data with complex experimental designs**. We extended the model to incorporate a flexible regularisation that enables the joint analysis of multiple omics as well as multiple sample groups (batches and/or experimental conditions). In addition, We implemented a GPU-accelerated stochastic variational inference framework, thus enabling the scalable analysis of potentially millions of samples. [Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02015-1) and [webpage](https://biofam.github.io/MOFA2/)

## License

http://creativecommons.org/licenses/by/4.0/

## Contact

ricard[dot]argelaguet[at]gmail.com

## Acknowledgements 

This thesis was written using the [template](https://github.com/kks32/phd-thesis-template ) created by Krishna Kumar