# Taxus-transcriptomic-study
Taxus, a medicinal plant, is widely recognized for its ability to natively produce paclitaxel. Placlitaxel, a plant sepcialized secondary metabolite, possesses tremendous anticancer functionality by stabilizing microtubules, thereby leading to death of cancer cells. Various types of cancers can be treated with paclitaxel, including ovarian and breast cancer. Paclitaxel has received its FDA approval in 1992. Since then, its market demand is rising.

Initailly, partial chemical syntheis has been employed for industrial production from its precursor extracted from the bark of a Taxus tree. However, due to the increased market demand of the compound, this production method requires excessive extraction of its precursor from the tree, posing serious threat to the survival of Taxus species. With the advancement of plant cell culture technology (PCC), it is possible to produce paclitaxel by developing Taxus PCC without desctroying its natural habitat. However, the continuous production of paclitaxel using Tasus PCC is limited by low yield. To tackle this issue, knowledge on global gene regulatory mechanisms responsible enhanced paclitaxel productio is necessary.

Here, we have conducted an extensive Taxus transcriptomic (RNA-seq) study on varying experimental condition to dicipher global gene regulatory mehcansim and discover critical genes might be responsible for enhanced paclitaxel production. This study will provide sufficient information for potential downstream engineering applications such as transforming Taxus cells for gene-overexpression, gene-knockout or RNAi, thereby leading to the generation of genetically modified cell-line that is cabaple of producing increased level of paclitaxel.

# Experimental design
From regular Taxus PCC, we have generated two phentoypes of different average aggregate size:
1. Small aggregates - Average diameter 500 $\mu m$
2. Big aggregates - Average diameter 1000 $\mu m$

We have then elicit the culture with methyl jasmonate (MeJA) to a final concentration of 150 $\mu M$
We have extracted RNA samples from all the cultures at two different time period.

**Table abbreviation:** **S** = Small aggregates, **B** = Big aggregates, **plus** = Elicited with MeJA, **min** = Mock elicited
| Sample name | Aggregate size | Treatment | Samplling time |
|-------------|----------------|-----------|----------------|
| S_18_plus | Small | Elicited with MeJA | 18 hour |
| S_72_plus | Small | Elicited with MeJA | 72 hour |
| S_18_min |	Small | Mock elicited |	18 hour |
| S_72_min |	Small |	Mock elicited |	72 hour |
| B_18_plus |	Big |	Elicited with MeJA | 18 hour |
| B_72_plus |	Big	| Elicited with MeJA |	72 hour |
| B_18_min |	Big	| Mock elicited |	18 hour |
| B_72_min |	Big |	Mock elicited |	72 hour |

# RNA-seq pipeline
After the RNA extraction, we have used Illumina Hi-Seq platform for read sequence. We have used several command line tool to process and align paired-end reads. Here is the list of tools we have used in this project to process the raw reads:

| Tools | Functions |
|-------|-----------|
| FastQC | Quality assessments of raw reads |
| TrimGalore | Remove low quality reads |
| Hisat2 | Read alignments to genome |
| Featurecounts | Calculation of transcripts abundance |
| DESeq2 | Calculation of differential gene expression |








