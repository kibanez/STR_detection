# STR_detection
Short tandem repeats detection in WGS

#### Motivation


Expansions, microsatellites or tandem repeats (short or long), are tracts of repetitive DNA sequences containing motifs ranging from 2-6 bases. Microsatellites are one of the most abundant type of variation in the human genome, after SNPs and Indels. But in contrast with these variations, microsatellites are rarely included or considered in WGS studies, in large due to a lack of tools capable of analysing them.



There are a number of regions in the human genome consisting of repetitions of short unit sequence (commonly a trimer). Such repeat regions can expand to a size much larger than the read length and thereby cause a disease. Fragile X Syndrome, ALS, and Huntington's Disease are well known examples. Many pathogenic repeat expansions have repeats spanning hundreds to thousands of base pairs. Thus it has been assumed that short-read sequencing technologies may not be able to accurately detect repeat expansions of pathogenic lengths. But for instance, as an exception, the `CGA` repeat in the HTT gene in which individuals with at least 40 repeat units (120bp) are at risk of developing Huntington's disease. For this repeat expansion, 150bp reads (V2) would be sufficient to identify the smallest pathogenic repeat expansions.


#### Methodology

Expansions are not routinely analysed in WGS studies, mainly due to a lack of adequate tools. First, not all reads that align to an STR locus are informative. If a pair-end read partially encompasses an STR locus, it provides only a lower bound on the number of repeats. Only reads that fully encompass an expansion can be used for exact expansion allelotyping. Second, mainstream aligners, such as BWA, generally exhibit a trade-off between run time and tolerance to indels. Profiling expansions variations (even for an expansion of 3 repeats in a trinucleotide repetition) would require a cumbersome gapped alignment step and lengthy processing times. Third, PCR amplification of an STR locus can create stutter noise, in which the DNA amplicons show false repeat lengths due to successive slippage events of DNA polymerase during amplification. Since PCR amplification is a standard step in library preparation for WGS, an expansion-detector strategy should explicitly model and attempt to remove this noise to enhance accuracy.

A number of methods have been developed to genotype microsatellites (Gelfand et al., 2104; Gymrek et al., 2012; Highman et al, 2013). We have tested lobSTR, popSTR and Expansion Hunter (EH) for this purpose. EH is the one with the best performance from the point of view of computational/cpu performance and the general outcome (with the negative and positives samples we analysed so far)*[]:


###### Expansion Hunter

Expansion Hunter aims to estimate sizes of such repeats by performing a targeted search through a BAM/CRAM file for reads that span, flank, and are fully contained in each repeat. It detects pathogenic repeat expansions from PCR-free, WGS data.

It identifies reads that either

* span the repeat (SPANNING reads)
* include the repeat and the flanking sequence on one side of the repeat (FLANKING reads)
* are fully within the repeat, called "in-repeat" reads (IRRs)

For short repeats (less than the read length of the sequence data) --> it calculates the repeat length using SPANNING and FLANKING reads.

For longer repeats (longer than the read length) --> it identifies and counts the IRRs



#### Validation

We need to see how EH estimates the tandem repeats in the genomes, and so, samples with known repeats (from the experimental point of view) are required

###### Positive samples - Coriell

The Mike Eberle's group at Illumina (San Diego) measured the FMR1 repeat expansionin 8 samples obtained from the Coriell biorepository (https://catalog.coriell.org/). FMR1 is a pathogenic repeat expansion where diagnosis depends on the repeat size: Tremos Ataxia Syndrome corresponds to a repeat consisting of 55-200 repeats (165-600bp) and Fragile X Syndrome is diagnosed with 200 or more repeats (>600bp). The experimentally detemined size of the FMR1 repeat expansion in these samples ranged from 76-645 repeats (228 - 1,935 bp).

They identified 20 off-target regions where `CGG` IRR pairs may disalign. They show in the paper that EH correctly categorised each sample as either pre-mutated (55-200 repeats) or expanded (>200 repeats), though it underestimated the repeat size in 3 of these samples.

We asked them for these samples and we got the corresponding BAMlets (mini BAM files) containing only the FMR1 gene (the one that already is 'affected'). We ran EH on these 8 samples and obtained the same outcome as the Illumina group.

###### Negative samples

We analysed the genomes coming from the NHNN (Henry Houlden and Arianna Tucci), that have been already analysed experimentally. We computed the performnce of EH for them, analysing SCA1,2,3,6 and 7.

