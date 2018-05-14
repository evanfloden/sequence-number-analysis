# embeded-analysis-nf
### An analysis to determine the effect of the inclusion of more sequences in multiple sequence alignments.


## Summary
This experiment illustrates the effect of sequence number on the accuracy of multiple sequence alignments.
It has been often assummed that the addition of more homologous sequences to a dataset results in 
alignments with higher accuracy owing to an increase in the available information content. Here we show
that after increasing the number of sequences beyond approximatly 1,000, the quality of alignment degrades.
We also show we can patially overcome this degredation using the regressive multiple sequence alignment procedure
where 'buckets' of sequences are aligned. 


## Replicating this Analysis

This analysis can be replicated in full. 

An interactive notebook writen in R Markdown provides complete instructions for doing so.

See the folder 'analysis' the base directory of this reposistory.


## Datasets

The datasets used are the original HomFam data. This dataset contains 94 families of protein. 
Each family has a set of reference sequeces which originate from structure in the PDB and a set
of sequences ranging from  94 sequences for the smallest dataset (seatoxin) up to 93,681 sequences 
(Retroviral aspartyl protease - rvp).


### Sequence Sets
For each familiy, sequence sets where created of the following sizes (up to the maximum number of sequences in that family): 

0 - 5 - 10 - 25 - 50 - 100 - 200 - 400 - 600 - 800 - 1,000 - 
1,250 - 1,500 -1,750 - 2,000- 2,500 - 3,000 - 4,000 - 5,000 - 
7,500 - 10,000 - 15,000 - 20,000 - 50,000 - 100,000

Each sequence set was generated 10 times.
For example, the sequence set with the filename `rvp.7500.1.fa` is the rvp familiy with 7,500 sequences, replicate 1.

### Guide Trees
The guide trees generated can be found in the `data/guide_trees` directory.

For example, the guide tree with filename `CLUSTALO.7500.1.dnd` is the guide tree generated with Clustal Omega for the `rvp.7500.1` sequence set.


### Alignments
The alignments generated can be found in the `data/alignments` directory.

There is three different types of alignments generated, `default`, `standard` and `regressive`.


* Default is the alignment generated when providing the MSA program wih the sequence set file (fasta) only.
* Standard is the alignment generated when providing the MSA program with the sequence set and a predetermined guide-tree.
* Regressive is the alignment generated using T-Coffee regressive framework upon providing the sequence set, a guide-tree and a bucket size.

For example, the alignment with the filename `rvp.7500.1.regressive.1000.CLUSTALO.with.CLUSTALO.tree.aln` is the alignment generated for the rvp dataset, with 7,500 sequences, replicate 1, using the regressive alignment framework with a bucket size of 1,000 and the Clustal Omega alignment method using the guide-tree from Clustal Omega.  

 
## Results

You can find the interactive notebook containing steps and results [here](rpubs link).

The non-interactive figure can be found in the figures directory.

Sequences, guide-trees and alignments can be found the 'data' directory of this repository.




