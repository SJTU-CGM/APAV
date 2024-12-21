
# APAV: An advanced pan-genome analysis and visualization toolkit for genomic pres-ence/absence variations
APAV is an advanced toolkit for comprehensive PAV analysis and visualization. It optimized the PAV detection process by incorporating an element-level analysis and adding PAV analysis for arbitrary regions in a genome. The resulted PAV profile can be viewed and checked in interactive web reports, providing researchers the ability to verify the read depth, region coverage and absent intervals of each PAV conveniently. Additionally, APAV offers various subsequent analysis and visualization functions based on the PAV table, including basic statistics, sample clustering, genome size estimation, and phenotype association analysis.


## Requirements
### PAV calling
* <b>Perl</b>
* <b>Samtools</b><br>Samtools is used to compute the read depth of positions in target regions. Please install it first and make sure it is under your `PATH`.

### Subsequent analysis and visualization based on the PAV profile
* <b>R</b> <br>R is utilized for visualization and statistical tests in the APAV toolbox. Please install R first and make sure `R` and `Rscript` are under your `PATH`.
* <b>APAVplot</b> <br>`APAVplot` is an R package specifically designed for visualization of PAV analysis. Follow the installation step, or you can install it by yourself. For the code and more details please see [here](https://github.com/SJTU-CGM/APAVplot).

## Installation procedures

1. You can download the APAV toolbox from github:
```
$ git clone https://github.com/SJTU-CGM/APAV.git
```
 Alternatively, you also could obtain the toolbox on the [APAV](https://cgm.sjtu.edu.cn/APAV/install.html) website and uncompress the APAV toolbox package:
```
$ tar zxvf APAV-v**.tar.gz
```

2. You need to add `apav` to `PATH` and add `lib/` to `PERL5LIB`
```
$ export PATH=$PATH:/path/to/APAV/:
$ export PERL5LIB=$PERL5LIB:/path/to/APAV/lib/:
```

3. The R package `APAVplot` is required for plotting
```
## Install "ComplexHeatmap" from BiocManager 
$ conda install r-BiocManager  ## OR: Rscript -e "install.packages('BiocManager')"
$ Rscript -e "BiocManager::install('ComplexHeatmap')" 
## OR install it from bioconda 
$ conda install bioconda::bioconductor-complexheatmap

## Install "APAVplot"
$ conda install r-devtools ## OR: Rscript -e "install.packages('devtools')" 
$ Rscript -e "devtools::install_github('SJTU-CGM/APAVplot')"
```
_You can skip this step if you are going to use APAVplot later to draw plots in the R environment_

4. Finally, you can test if the APAV toolkit is installed successfully by:
```
$ apav
```

If you could see the following content, congratulations! APAV toolkit is successfully installed. If not, see if all the requirements are satisfied; or you may contact the authors for help.
```
Usage: apav <command> ...
Available commands:

    Pipeline:
        geneBatch               Automatically execute main commands for genes
        generalBatch            Automatically execute main commands for the general target regions

    Extract positions:
        gff2bed                 Extract the coordinates of target regions from a GFF format file
        
    Calculate coverage:
        staCov                  Calculate coverage of target regions
        mergeElecov             Merge neighboring elements with the same coverage
        covPlotHeat             Plot a heatmap to give an overview of the coverage profile across samples
        
    Determine PAV:
        callPAV                 Determine presence/absence variations based on coverage
        gFamPAV                 Determine gene family presence/absence based on the gene PAV table
        mergeElePAV             Merge neighboring elements with the same PAV

    Estimate genome size:
        pavSize			Estimate the size of pan-genome and core-genome from the PAV table
        pavPlotSize		Draw estimated growth curves

    PAV analysis:
        pavPlotStat             Plot a half-violin chart to show the number of regions in each group of samples
        pavPlotHist             Plot a ring chart and a histogram to show the classifications and distribution of target regions
        pavPlotHeat             Plot a complex heat map to give an overview of the PAV profile
        pavPlotBar              Plot a stacked bar chart to show the classifications of target regions in all samples
        pavPCA                  Perform PCA analysis for the PAV table and plot results
        pavCluster              Cluster samples based on the PAV table and plot results

    Phenotype assocation:
        pavStaPheno             Perform Fisher's exact test and Wilcoxon tests to determine phenotype association
        pavPlotPhenoHeat        Show the main result of phenotype association analysis with a heat map
        pavPlotPhenoBlock       Display the percentage of samples containing target regions in each group of a discrete phenotype
        pavPlotPhenoMan         Draw a Manhattan plot to show the results of a given phenotype
        pavPlotPhenoBar         Show the relationship between a specific genomic region and a specific phenotype in a bar plot
        pavPlotPhenoVio         Show the relationship between a specific genomic region and a specific phenotype in a violin plot

     Visualization of element regions:
        elePlotCov              Display the coverage of elements in a specific target region
        elePlotPAV              Display the PAV of elements in a specific target region
        elePlotDepth            Display the depth of elements in a specific target region
```

The usage information for each command can be shown with the `--help` or `-h` option after each command name. You also could get parameter lists on the [APAV](https://cgm.sjtu.edu.cn/APAV/usage.html) website.

## Quick start

### Integrated commands

The main steps can be automatically executed with the `geneBatch`/`generalBatch` command. The `geneBatch` command is used for gene region and takes the GFF file as input. The `generalBatch` command is used for any target region and takes the BED file as input.

```
cd ${APAV_PATH}/demo/

## demo1: Some genes on human chromosome 19
apav geneBatch --gff demo1_gene.gff3 --bamdir bam --pheno demo_sample.pheno --fa demo.fa.gz --fam demo1_gene.fam --up_n 10 --down_n 10 --chrl demo1.chrl --out demo1

## demo2: Some proteins on human chromosome 19
apav generalBatch --bed demo2_general.bed --bamdir bam --pheno demo_sample.pheno --fa demo.fa.gz --out demo2

## demo3: Some repeat sequences on human chromosome 19
apav generalBatch --bed demo3_general.bed --bamdir bam --pheno demo_sample.pheno --fa demo.fa.gz --rmele --out demo3
```

The PAV reports for demos can also be viewed [here](https://cgm.sjtu.edu.cn/APAV/demo.html).

### Step-by-step commands

#### 1. Coordinate extraction
For genes, you need to use the `gff2bed` command to extract the coordinates of genes, genetic elements and bins in upstream/downstream. It merges elements with the same coordinates and outputs them in BED format. For general target region, skip this step.
```
apav gff2bed --gff demo1_gene.gff3 --out demo1.bed
```
Options `--chrl`, `--up_n`, `--up_bin`, `--down_n`, `--dwon_bin` are used for extracting upstream and downstream elements.
```
apav gff2bed --gff demo1_gene.gff3 --chrl demo1.chrl --up_n 10 --up_bin 100 --down_n 10 --down_bin 100 --out demo1.bed
```			

#### 2. Coverage calculation
Use the `staCov` command to compute the coverage of regions. It counts the percentage of covered bases in the whole region and each element region.
```
apav staCov --bed demo1.bed --bamdir bam --asgene
```
The `covPlotHeat` command will plot a heatmap to give an overview of coverage profile.
```
apav covPlotHeat --cov demo1.cov
```
		
#### 3. PAV determination
Based on the coverage, the `callPAV` command determines the presence-absence variation. It generates PAV profiles and two interactive web reports.
```
apav callPAV --cov demo1.cov --pheno demo_sample.pheno --fa demo.fa.gz --gff demo1_gene.gff3 
apav callPAV --cov demo1_ele.cov --pheno demo_sample.pheno
```

#### 4. Gene family PAV determination
For genes, the `gFamPAV` command allows further determination of gene family PAV profile.  
```
apav gFamPAV --pav demo1_all.pav --fam demo1_gene.fam
```
	
#### 5. Genome size estimation
Based on the PAV table, you can use the `pavSize` command to estimate genome size by simulating the size of the pan-genome and core-genome.
```
apav pavSize --pav demo1_all.pav
```
```
## Estimation in groups
cat demo_sample.pheno | cut -f 1,2 > demo_sample.group
apav pavSize --pav demo1_all.pav --group demo_sample.group
```
The `pavPlotSize` command can draw the growth curve of genome estimate.
```
apav pavPlotSize --size demo1_all.size
```

#### 6. Common PAV analysis and visulization
APAV provides various commands for common PAV analysis. The `pavPlotStat` command shows the total number of regions in all samples. The `pavPlotHist` command shows the classifications and distribution of regions. The `pavPlotHeat` command gives an overview of the PAV table. The `pavPlotBar` command shows the composition of each sample. The `pavPCA` command performs PCA analysis. The `pavCluster` command clusters samples based on the PAV table.
```
apav pavPlotStat --pav demo1_all.pav
apav pavPlotHist --pav demo1_all.pav
apav pavPlotHeat --pav demo1_all.pav
apav pavPlotBar --pav demo1_all.pav
apav pavPCA --pav demo1_all.pav
apav pavCluster --pav demo1_all.pav
```

#### 7. Phenotype association analysis
Use the `pavStaPheno` command to determine phenotype association.
``` 
apav pavStaPheno --pav demo1_all.pav --pheno demo_sample.pheno
```
The `pavStaPhenoHeat` command gives an overview of significantly phenotype-related regions. The `pavPlotPhenoBlock` command is used to display discrete phenotype. The `pavPlotPhenoMan` command draws a Manhattan plot. The `pavPlotPhenoBar` and `pavPlotPhenoVio` commands show the relationship between a certain genomic region and a certain phenotype.
```
apav pavPlotPhenoHeat --pav demo1_all.pav --pheno_res demo1_all.phenores
apav pavPlotPhenoBlock --pav demo1_all.pav --pheno demo_sample.pheno --pheno_res demo1_all.phenores --pheno_name Gender
apav pavPlotPhenoMan --pav demo1_all.pav --pheno demo_sample.pheno --pheno_res demo1_all.phenores --pheno_name Gender
apav pavPlotPhenoBar --pav demo1_all.pav --pheno demo_sample.pheno --pheno_name Location --region_name ENSG00000233493.3
apav pavPlotPhenoVio --pav demo1_all.pav --pheno demo_sample.pheno --pheno_name Age --region_name ENSG00000254415.3
```
These steps also apply to elements.
	
#### 8. Visualization of element regions
For the focused target region, you can use the `elePlotCov`/`elePlotPAV` command to observe the coverage/PAV of elements. Furthermore, the `elePlotDepth` command can display the read depth in target regions.
```
grep 'ENSG00000126251.6' demo1_gene.gff3 > ENSG00000126251.6.gff3
grep -E 'Annotation|ENSG00000126251.6' demo1_ele.cov > ENSG00000126251.6.elecov
grep -E 'Annotation|ENSG00000126251.6' demo1_ele_all.pav > ENSG00000126251.6.elepav
apav elePlotCov --elecov ENSG00000126251.6.elecov --pheno demo_sample.pheno --gff ENSG00000126251.6.gff3
apav elePlotPAV --elepav ENSG00000126251.6.elepav --pheno demo_sample.pheno --gff ENSG00000126251.6.gff3
apav elePlotDepth --ele ENSG00000126251.6.elecov  --bamdir bam --pheno demo_sample.pheno --gff ENSG00000126251.6.gff3
```
