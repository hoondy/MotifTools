# MotifTools
A collection of python tools to evaluate a variant D-score (motif-breaking or motif-gaining power) on TF motif

![MotifTools](https://docs.google.com/a/yale.edu/uc?authuser=1&id=0B1hrcjjDSLuXamFwakx1SHpuc28&export=download)

## Dependencies

* Python 2.6 or higher (tested on Python 2.7.12)
* BioPython 1.66 or higher
* NumPy
* BedTools

## Config

1. Make a copy of config.ini from config-sample.ini
2. Specify path to bedtools
3. Specify path to JASPAR database pfm file (i.e., JASPAR TF profiles, 2016 core non-redundant vertebrates, can be downloaded [here](http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_vertebrates.txt))
4. C_PRECISION tells how many digits to preserve when calculating P-value
5. C_PSEUDOCOUNTS is added to PFM to avoid division by zero
6. C_PVAL_THRESHOLD is P-value threshold for reporting TF motif significance on either reference or alternate genome. For example, suppose C_PVAL_THRESHOLD is 0.0001 and there is a variant that causes TF motif P-value to change from 0.001 to 0.1. Although it has D-score of 20, it does not get reported because neither TF motif on reference genome nor on alternate genome is significant. If a TF motif has P-value of 0.001 on reference genome and if a variant makes a TF motif to have P-value of 0.0001, it will be reported with D-score of -10.

## Inputs

* Peak BED file
* Variant VCF file
* Reference genome FASTA file
  * For example, hg19 genome can be downloaded from [here](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/)  

## Usage

```shell
python runDscore.py --sample {"MCF-7"} --tf {"CTCF"} --bed {peak.bed} --vcf {variant.vcf} --ref {reference_genome.fa}
```

or

```shell
python runDscore.py -s {"MCF-7"} -t {"CTCF"} -b {peak.bed} -v {variant.vcf} -r {reference_genome.fa}
```

## Output

Output is saved in BED format (i.e., genomic coordinate is 0-based).

```
chr1	205304887	205304906	chr1:205304888-205304906_+	12.5378302556	+	3.20285835187e-06	5.74540645175e-05	TCCTCA[C]TAGGGGGCAGCA	TCCTCA[T]TAGGGGGCAGCA	15.122657723	9.53769650652	7
chr14	107239045	107239065	chr14:107239046-107239064_-	-0.857778555555	-	7.40751347621e-05	6.07987421972e-05	AGGACATGAGGTGGC[A]CTC	AGGACATGAGGTGGC[G]CTC	8.96464812219	9.41133123776	16
chr1	115499385	115499405	chr1:115499386-115499404_-	1.22679468983	-	6.97545074217e-05	9.25234344322e-05	tcagcactagaggatg[C]tc	tcagcactagaggatg[G]tc	9.10192102726	8.45112655346	17
chr15	79888612	79888632	chr15:79888613-79888631_-	-1.33949922329	-	8.03446164355e-06	5.90210402152e-06	CAGCCTCTAGGAGGAGCT[C]	CAGCCTCTAGGAGGAGCT[G]	13.5124960611	14.0679325841	19
chr20	57023063	57023083	chr20:57023064-57023082_-	-1.16307241472	-	6.53455936117e-05	4.99929847138e-05	CAGACAGCCGGGGGCAGG[G]	CAGACAGCCGGGGGCAGG[A]	9.24950881325	9.84527805058	19
chr1	201528952	201528971	chr1:201528953-201528971_+	12.2863623638	+	7.54283973947e-07	1.27694183902e-05	CCGCCACC[A]GAGGGTGCAG	CCGCCACC[G]GAGGGTGCAG	17.3776446873	12.6446297995	9
chr1	201528952	201528971	chr1:201528953-201528971_+	-4.59613738195	+	7.54283973947e-07	2.61770765064e-07	CCGCCACCAGAGGGTGC[A]G	CCGCCACCAGAGGGTGC[C]G	17.3776446873	18.8381967673	18
chr5	176541052	176541071	chr5:176541053-176541071_+	1.43157045253	+	8.56964798004e-05	0.000119157128211	TAG[T]CACACGGTGGCGACA	TAG[G]CACACGGTGGCGACA	8.62883752903	7.85123023661	4
chr17	55825082	55825102	chr17:55825083-55825101_-	-7.89814065592	-	8.06438401924e-06	1.30845000967e-06	AGACCAGTAGGTGTCA[T]CA	AGACCAGTAGGTGTCA[G]CA	13.5051877266	16.5578288634	17
chr17	53591549	53591569	chr17:53591550-53591568_-	3.72653409279	-	1.96814653464e-09	4.64206095785e-09	TGGCCACTAGATGGCACC[A]	TGGCCACTAGATGGCACC[G]	23.6278904949	23.0321212576	19
```

* COLUMN 01: chromosome
* COLUMN 02: start position of motif
* COLUMN 03: end position of motif
* COLUMN 04: 1-based genomic coordinate ID
* COLUMN 05: D-score
* COLUMN 06: strand
* COLUMN 07: P-value on reference genome
* COLUMN 08: P-value on alternate genome
* COLUMN 09: motif sequence on reference genome
* COLUMN 10: motif sequence on alternate genome
* COLUMN 11: -10 * log(P-val_Ref)
* COLUMN 12: -10 * log(P-val_Alt)
* COLUMN 13: variant position w.r.t. motif

## How to interpret D-score

D-score is a motif "[D]isruptive score" of a variant. It is calculated by difference between P-value between reference genome and alternate genome.

```
D-score = [-10 * log(P-val_Ref)] - [-10 * log(P-val_Alt)]
D-score = -10 * log(P-val_Ref/P-val_Alt)
```

* *Positive D-score denotes a variant is decreasing the likelihood of TF to bind the motif (motif-break)*
* *Negative D-score denotes a variant is increasing the likelihood of TF to bind the motif (motif-gain)*


## Future Improvements

* Current version has been optimized for calculating D-score with equal nucleotide background (A:C:G:T=1:1:1:1). While changing the proportion of nucleotide background is supported in current version, it can be slow and precision maybe lost. Our internal test has indicated that changing the proportion of nucleotide background has very minimal impact on P-value. The future version should address this issue.
* Motif logo of BEFORE and AFTER
