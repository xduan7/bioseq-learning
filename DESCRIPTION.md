# Biological Sequence Learning with Conserved Domain


## Table of Contents
-   [Motivation and Rationale](#motivation-and-rationale)
-   [Input Format and Scale](#input-format-and-scale)
-   [Sequence Data Processing](#sequence-data-processing)
-   [Predictive Target](#predictive-target)
-   [Learning Model](#learning-model)
-   [Experiment Design](#experiment-design)
-   [Results and Analysis](#results-and-analysis)
-   [Conclusion](#conclusion)


---
## Motivation and Rationale 


### What is conserved domain and CDD

[Wikipedia: Conserved Domain Database](https://en.wikipedia.org/wiki/Conserved_Domain_Database):
> Domains can be thought of as distinct functional and/or structural units of a protein. These two classifications coincide rather often, as a matter of fact, and what is found as an independently folding unit of a polypeptide chain also carries specific function. Domains are often identified as recurring (sequence or structure) units, which may exist in various contexts. In molecular evolution such domains may have been utilized as building blocks, and may have been recombined in different arrangements to modulate protein function. CDD defines conserved domains as recurring units in molecular evolution, the extents of which can be determined by sequence and structure analysis.
> The goal of the NCBI conserved domain curation project is to provide database users with insights into how patterns of residue conservation and divergence in a family relate to functional properties, and to provide useful links to more detailed information that may help to understand those sequence/structure/function relationships. To do this, CDD Curators include the following types of information in order to supplement and enrich the traditional multiple sequence alignments that form the foundation of domain models: 3-dimensional structures and conserved core motifs, conserved features/sites, phylogenetic organization, links to electronic literature resources.


### Why Conserved Domain as Prediction Target
-   [ ] what can this model do other than predict CD
-   [ ] what does this model understand
-   [ ] any other long-term goals


---
## Input Format and Scale


### Nucleotide or Protein Sequences as Input 
-   [ ] advantages and disadvantages of input choices


### Type and Scale of Sequences for Training
-   [ ] on the smaller scale
-   [ ] on the bigger scale 
-   [ ] ultimately ... 


---
## Sequence Data Processing

### CD Search 
-   [ ] which CD database(s) to use? There are around 16,000 NCBI-curated CDs, while other 40,000 CDs are sourced from elsewhere. Use NCBI-curated CDs only will limit the search, which means less computation during the training and less knowledgable for trained the model
-   [ ] computation feasibility

### Other Input Annotations and the Potential Usage
-   [ ] annotation for NCBI-curated CDs. There are binding sites, loop/helix and other region annotation ...
-   [ ] protein family label
-   [ ] other annotations


---
## Predictive Target

-  [ ] choose between the following predictive targets for training

* predict the start/center/end of CDs - translation
* predict the belonged CDs(s) for each and every monomer (N-to-N) probably onluy suitable for small number of PSSMs ... and short sequence
    * training as a multi-label task, then demonstrating the distribution of the latent space of CDs, showing that some of the similar CDs are much closer
    * training as a multi-label task, but the errors are weighted with similarities between the predicted CDs and the true ones
    * training with hierarchical information of CDs (only viable for NCBI curated CDs)
* predict the belonged CD(s) for a specific word (which CD it's in and where it's at relatively)
* CD as input and only predict word and the CD it's in
* predict the whole masked CD with given genome context (or the next CD)


Input choice: gene or genome
Output choice: masked monomer, masked CD, masked segment, or CD sequence



---
## Learning Model

-   [ ] choose between the following NLP models

* BERT
* Transformer-XL
* GPT/minGPT


---
## Experiment Design


---
## Results and Analysis


---
## Conclusion

