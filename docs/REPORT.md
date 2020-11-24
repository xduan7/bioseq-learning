# Biological Sequence Learning with Deep Learning


## Table of Contents
-   [Introduction](#introduction)
-   [Data Hierarchy](#data-hierarchy)
    -   [Biological Segmentation](#biological-segmentation)
    -   [Statistical Segmentation](#statistical-segmentation)
-   [Learning Model](#learning-model)
    -   [Network Architecture](#network-architecture)
    -   [Training Task](#training-task)
-   [Evaluation Method](#evaluation-method)


---
## Introduction

This project applies the state-of-the-art deep learning models to biological sequences (amino acid or nucleotide sequences). 
There are several notable challenges for learning on biological sequences:
-   **biological sequences, especially nucleotide ones, could be extremely long**
-   **knowledge is buried deep inside the sequences, and therefore hard to learn**
-   the sheer amount of data and the quality issues;

The first two challenges are more relevant to the research of machine learning and its application in biology, and interestingly enough, despite being two distinct challenges, they could be solved at the same time with [a better model](#learning-model) and/or [data hierarchy](#data-hierarchy).


---
## Data Hierarchy

In the context of this project, the hierarchy for biological sequences indicates the inner relation and structure of smaller segments. For instance, the protein coding regions on top of genome sequences is a level of hierarchy.

The hierarchical structure of the biological sequences makes the learning easier in two ways: 
1.  sequences are presented in a more structured and hopefully more "understandable" way for machine learning models
2.  longer sequences could be segmented into shorter chunks, which are easier to process and learn

The segmentation could be either [biological](#biological-segmentation) or [statistical](#statistical-segmentation), as explained separably in the subsections below.


### Biological Segmentation

This segmentation method is based on the human understanding of biological sequences. 
For example, a genome can be segmented into coding regions and non-coding regions for proteins; and on coding regions, there could be [conserved domains]((https://en.wikipedia.org/wiki/Conserved_Domain_Database)), which is an extra level of hierarchy. 
Another good example is the [secondary structure](https://en.wikipedia.org/wiki/Protein_secondary_structure) of protein, which segments the sequences into alpha-helix, beta-sheet, etc. 

The applied biological knowledge for this segmentation method is very likely to benefit the learning process. However, the segmented sequence could still be too long for learning. Some coding regions in E. Coli genomes have more than 5,000 nucleotide bases, which cannot easily fit into some deep learning models, especially the ones that rely on self-attention.

#### Conserved Domain as Biological Segmentation

Conserved domains "can be thought of as distinct functional and/or structural units of a protein". Ideally, if one can segment the genomes or proteins sequences based on conserved domain, the data hierarchy might provide the model with valuable information. 

However, the "distribution", for lack of a better word, of different conserved domains are rather skewed. There are more than 13,000 different conserved domains found in E. Coli genomes, but most of them are rather rare, appearing less then 500 times across more than 1,000 whole genomes. So training with conserved domains as if they are words in natural languages is probably not feasible. 
Not to mention that the occurrences of conserved domains are relatively rare already. 
And the definition of conserved domain is not clear-cut, but largely defined on the database to perform conserved domain search with.

Combine all the factors together, conserved domains are probably not the best way to provide hierarchy. If we are going to try segmenting the sequences with conserved domains anyway, we should focus only on the most common conserved domains from our dataset.

### Statistical Segmentation

Unlike biological segmentation, this segmentation method is purely based on the statistical property of the sequences. 
For example, one can segment a long genome sequence into multiple smaller ones within some pre-defined length. Another good example is to segment a genome sequence with the occurrence of start and stop codon.

The biggest disadvantage of this method is that it provides little if none biological knowledge to the learning model, but it is easily to implement. Moreover, the segmented sequences could be shorter and faster to learn.


--- 
## Learning Model

A model suitable for the learning task of biological sequences should be able to digest ultra-long sequences and provided good performance at the same time. The definition of learning model can be boiled down to two orthogonal components: the [network architecture](#network-architecture) and the [training task](#training-task).

### Network Architecture

There are several network architecture suitable for learning on biological sequences:
-   dense neural network (baseline)
-   **1-D convolution neural network**
-   transformer (BERT) encoder
    - **transformer (BERT) encoder with sparse attention**

Out of all the models listed above, dense models and the original transformer (with dense self-attention) are not as suitable as the other two because of the length of the biological sequences. Some level of sparsity like convolution might be the only solution.
Note that recurrent neural network might not be a good idea here because the information might be lost as the recurrent process continues on the ultra-long sequences.


### Training Task

Since there is an enormous amount of genome/protein sequences available without too much informative labeling or easy supervised-learning task, it is only reasonable to apply unsupervised learning.

**TODO: research on the unsupervised tasks for sequences.**


--- 
## Evaluation Method
<!---  
-   [ ] tSNE visualization of genetic codes
-   [ ] faster conserved domain search (BLAST score)
-   [ ] protein family/function classification
-   [ ] protein structural similarity (TM alignment score) prediction 
-   [ ] mutation effect
-   [ ] gene regulatory network prediction inference
-   [ ] co-varying mutations 
--->

