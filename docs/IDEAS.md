This file contains ideas and thoughts that are not fully baked but worth going back and examine. 


## Hierarchical Learning for Genome Sequences
There are four levels of learning that we could perform 

###  character level (ATGCs or AAs) sequence -> protein domains:

This is more or less like an strictly ordered translation; or some character-level language to words with semantic meaning

- how to decide on the input sequence length (start and end); could it be the same as the automated gene identification process?
- are the input sequences made of nucleotide bases or amino acids? if it's the latter, then are we not using non-coding regions on genomes (since they are not going to be translated into AAs and there is no reading frame)?
- will the non-domain regions be somehow coded into the output sequence (from a symbol to a learned vector)? or do we just ditch them all together?

Note that even if this step failed, we are still going to have a training set of protein domain sequences for the next step, which are extracted from existing software. The reasons for this step are (1) same protein domains can have different vector representation when there are some sequencing differences, and (2) regions that are not domains COULD have a representation of some sort, but it's all up to implementation.

Question: are there any other units for protein/gene sequences? like protein subunit, etc.

### protein domains -> protein functions

### masked protein (domain) prediction on the whole genome

### genotype -> phenotype