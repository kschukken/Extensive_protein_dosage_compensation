# Extensive_protein_dosage_compensation
Updated repository for code and intermediate data files for "Extensive dosage compensation in aneuploid human cancers". Update includes ovarian tumor analysis and additional control datasets. 

Depository of all R files for “Extensive protein dosage compensation in aneuploid human cancers”
This is a repository for the BioRxiv preprint "Extensive protein dosage compensation in aneuploid human cancers". By Klaske M Schukken and Jason M Sheltzer. This is the updated repository, with additional files that analyze Ovarian tumor data and additional control datasets. 
The Summary_of_files.xlsx is an overview of which R file generates which datasets and graphs, also listed below.
Many of the code files need additional datasets retrieved from outside sources-- The source of the datasets are often listed in the ## Step 1: Get data ## section of the R files. If data source is not listed in the R code itself, it will be credited in the "Extensive protein dosage compensation in aneuploid human cancers" BioRXiv pre-print. Several intermediate datasets are made available for ease of use.
All code was build on R version 3.6.3 GUI 1.70 El Capitan build, and R studio version 1.1.383.
In order to run the code, open the .R files, download the indicated datasets, set the proper working directory, and run through the code. Some code may take several days to run. Often, datasets that take a long time to generate are also available as an intermediate file. 


ORGANIZATION
R CODE:
Protein_RNA_filtered_CellLine.R
·      All figures
·      Filter protein expression and RNA expression data to isolate cell lines with RNA expression, protein expression, and arm call data
·      Generates (among other files): Protein_Expression_filtered.csv

Protein_RNA_expression.PerCell_v2.R
·      Figure 2-6, S2, S4-7
·      Difference in gene expression when that gene is located on an aneuploid chromosome (gain/loss relative to neutral copy number). 
·      Also plots RNA or protein expression by chromosome copy number category (gain, neutral, loss) boxplots
·      Generates (among other files): CN.Diff.xRNA.yProt.ThreeGroups

Protein_filtered_analysis_3cat.R
·      Figure 1, S1
·      Calculate protein expression difference for all genes (not just genes on aneuploid chromosome) upon chrm arm gain or loss (ex. gain chrm 5, protein expression diff in chrm 1,2,3,4,5,6, etc. )
·      Generates (among other files): ChrmArmGain.5q.Protein.Diff.all_min10cells.csv

RNA_filtered_analysis_3cat.R
·      Figure 1, S1
·      Calculates RNA expression difference for all genes (not just genes on aneuploid chromosome) upon chrm arm gain or loss (ex. gain chrm 5, RNA expression diff in chrm 1,2,3,4,5,6, etc. )
·      Generates (among other files): ChrmArmGain.5q.RNA.Diff.all_min10cells.csv

ChrmArm.Scatterplot.Protein.RNA.R
·      Figure 1
·      Protein and RNA expression difference per chrm arm scatterplot (ex. Gain chromosome 5q: protein expression difference on chrm 5q, and on other chromosomes)

DNA_Protein_RNA_filtered.R
·      Figure 1, S1
·      Make graphs of chrm arm difference from neutral-ploidy, for specific cells (ex. MCF7)

Protein_RNA_Heatmap_filtered_min10percategory.R
·      Figure 1
·      Generate heatmaps: protein/RNA expression difference upon chrm arm gain/loss, and additional summary heatmap

Protein_RNA_Corr_v2.R
·      Figure 2, 6
·      Graph: Protein or RNA expression per cell in gain/neutral/loss categories

Protein_RNA_1d.plot_v2
·      Figure 3, 6, S4, S7
·      Plot 1dimentional density plot of gene expression changes per group (Protein and RNA). Also look at low vs high aneuploidy cells only

gProfile_Quantile_Bargraph_v2.R
·      Figure 3, 4, S4 ,S6, S7
·      Make graphs of key Gene Ontology enrichment terms (for multiple analysis')

Yeast_proteomics_v2.R
·      Figure S5
·      Aneuploid yeast analysis (Protein & RNA)

DownSyndrome_analysis_v2.R
·      Figure S5
·      Down syndrome fibroblast analysis (Protein and RNA)

Stingele_Protein_filter_v2.R
·      Figure S5
·      Stable aneuploid cell line analysis (Protein and RNA)

RPPA.Difference.R
·      Figure S2, S3
·      RPPA protein expression data analysis

Aneuploid_Score.R
·      Figure 4, S6
·      Generate cellular aneuploidy scores
·      Generates (among other files): Aneuploidy.Score.perCell.xlsx

Protein_expression_data_GeneScore.R
·      Figure 4
·      Protein expression correlate with cellular aneuploidy, and high/low aneuploid cell only analysis
·      Generates (among other files): Protein_Corr_pvalue_GeneAneuploidy.csv

CCLE_RNA_Analysis_geneAn.R
·      Figure S6
·      RNA aneuploidy score correlation analysis
·      Generates (among other files): CCLE_RNA_Score2_GeneAneuploidy.csv

Variance_aneuploidy.R
·      Figure 5
·      Calculate RNA and Protein neutral-ploidy variance

Protein_buffering_factors_v2.R
·      Figure 5
·      Protein buffering factors. Calculate ROC AUC and make boxplots of scores/values
·      Generates (among other files): Protein.AllFactors.csv

TSG.OG_CNV_Difference_v2.R
·      Figure 6
·      Oncogene and tumor suppressor gene difference upon chromosome copy number changes and gene copy number changes

TCGA_ProteinDiff_Ovarian.R
·      Figure 7, S8
·      Ovarian tumor data analysis. Generating protein and RNA expression differences upon chromosome gain and loss, and plotting the corresponding data. 




ADDITIONAL FILES:

arm_calls_data_bendavid.csv
·      Aneuploid arm calls per chromosome arm per cell line analyzed. Used to generate aneuploidy score per cell line, and used to identify which cell lines have chromosome arm gains and losses per chromosome arm. Data from Cohen-Sharir et al. 2021 Nature, from the lab of Dr. Uri Ben-David.
 
bp_per_arm.csv
·      number of basepairs per chromosome arm. Used to generate a basepair-based aneuploidy score
 
CellLines_ProteinDosageCompensationManuscript.csv
·      List of cancer cell lines used for our CCLE data analysis. These cells have RNA expression, mass spec protein abundance and chromosome arm call data available. 

Score2.aneuploid.cell_line.csv
·      Cellular aneuploidy scores per cell line. Multupiple forms of the aneuploidy score are calculated, the gene_ploidy_score was used in the manuscript unless otherwise indicated.
 
RNA_Protein_Shared_cells.csv
·      Cell lines with both RNA expression data and protein expression data. Not all of these cells are also present in the Cohen-Sharir arm call data.
 
Protein_location_info.csv
·      Information about the chromosome arm location per gene.
 
Protein_ID_info.csv
·      Protein_IDs matched with gene symbol and uniprot accesion
 
RNA.Protein_loss.neutral.gain_Difference_Pvalue_min10points.csv
·      Dataset of gene expression difference on cells with chromosome gain or loss, relative to cells with a neutral ploidy for chromosome arm a gene is located on.
 
RNA.Protein_loss.neutral.gain_Difference_Pvalue_min10points_3cat.csv
·      Dataset of gene expression difference on cells with chromosome gain or loss, relative to cells with a neutral ploidy for chromosome arm a gene is located on.
with gene classifications as either "Anti-Scaling", "Buffering" or "Scaling" upon chromosome gain and loss
 
Trisomy 21 comparison.xlsx
·      Protein expression difference between Down Syndrome (Ts21) fibroblasts and control fibroblasts, set relative to cancer cell line protein expression difference. Yeast data from: Liu et al. 2017, Nature Communications. Systematic proteome and proteostasis profiling in human Trisomy 21 fibroblast cells
 
yeast-human aneuploidy.xlsx
·      Yeast protein expression difference upon chromosome gain. Data originated from: Dephoure et al. 2014, eLife; 3:e03023 Quantitative proteomic analysis reveals posttranslational responses to aneuploidy in yeast
 
Aneuploidy.Score.perCell.xlsx
·      Cellular aneuploidy Scores, calculated based on number of genes in aneuploid chromosomes. 

AmplificationRatio_perGene1.75.csv
·      Calculation of the percent of cells in the filtered dataset with gene amplifications per gene. Used in initial factor analysis, but not included in the final manuscript. Data from CCLE gene copy number data.
 
Protein.AllFactors.csv
·      CCLE data: A list of all genes, difference upon chromosome gain or loss, and factor scores/values per gene. 5' UTR length and 3' UTR length were not included, because there are frequently multiple 5' or 3' UTR lengths per gene, and this duplication can scew data when investigating other factors. See Protein_buffering_factors_v2.R, to incorperate 5' UTR length and 3'UTR length if desired.
 
Protein.AllFactors_Ovarian.csv
·      Ovarian tumor data: A list of all genes, with difference upon chromosome gain or loss, and factor scores/values per gene. 5' UTR length and 3' UTR length were not included, because there are frequently multiple 5' or 3' UTR lengths per gene, and this duplication can scew data when investigating other factors. See Protein_buffering_factors_v2.R, to incorporate 5' UTR length and 3'UTR length if desired.
 
Factor.ROCAUC.correlationscore.csv
·      A list of all investigated factors, their AUC ROC values with regard to buffering, scaling and anti-scaling genes, and the correlation coefficients with RNA and Protein expression difference.
 
TSG.ExpressionDifference.Gain.Loss.csv
·      Protein and RNA expression differences for tumor suppressor genes

Oncogene.ExpressionDifference.Gain.Loss.csv
·      Protein and RNA expression differences for Oncogenes

RPPA.MassSpec.difference.csv
·      Difference in protein expression when genes are located in aneuploid chromosome vs neutral chromosomes for data from mass spec and RPPA generated proteomics. 
 
TCGA_ovarian.diff.pvalue.RNA.Protein.ThreGroups.csv
·      Ovarian tumor data: Dataset of gene expression difference on cells with chromosome gain or loss, relative to cells with a neutral ploidy for chromosome arm a gene is located on. with gene classifications as either "Anti-Scaling", "Buffering" or "Scaling" upon chromosome gain and loss
 
Supplemental Data 7 gene lists
·      List of genes per go term used. You can also download GO term gene lists from the gene ontology database directly. 


Control datasets:
RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_highRNA_Min10Cells.csv
·      Control dataset looking at protein and RNA expression differences when all genes with the lowest 20% of RNA expression have been removed. 

RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_Lung.min10cells.csv
·      Control dataset looking at protein and RNA expression differences from only NSCLC cell lines. 

RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_noFlip_Min10Cells.csv
·      Control dataset looking at protein and RNA expression differences after removing all genes whole gene copy numbers are gained when their chromosome are lost, or whose gene copy number is decreased when their chromosome copy number is lost. 
 
RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_NoLowReplication_Min10Cells.csv
·      Control dataset looking at protein and RNA expression differences when all genes with the lowest 20% reliability scores have been removed. Also all genes without a reproducibility score were removed.  
 
RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_noMut_Min10Cells.csv
·      Control dataset looking at protein and RNA expression differences when all datapoints with mutations listed in the Depmap Mutation database were removed. 
 
RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_noMut_noFlip_noLowRNA_nolowRep_min10Cells.csv
·      Control dataset looking at protein and RNA expression differences when all mutated genes, “flipped” genes, low RNA expressing genes and low reproducibility genes have been removed.  The “Merged-control” dataset.  
