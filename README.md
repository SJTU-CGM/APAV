
# APAV: An advanced pan-genome analysis and visualization toolkit for genomic pres-ence/absence variations
APAV is an advanced pan-genome analysis and visualization toolkit for genomic presence-absence variations. It accepts the GFF file of genes and the BED file of any target regions. It takes map-to-pan strategy and computes coverage of target regions at the whole level and element level. You can review and check all PAVs and samples in the automatically generated interactive web reports. Based on the PAV table, APAV also offered various subsequent analysis and visualization functions, including basic statistics, sample clustering, genome size estimation, and phenotype association analysis.



## Requirements
### PAV calling
* Perl
* Samtools<br>Samtools is used to compute the read depth of positions in target regions. Please install it first and make sure it is under your PATH.

### Subsequent analysis and visualization based on the PAV profile
* R <br>R is utilized for visualization and statistical tests in the APAV toolbox. Please install R first and make sure R and Rscript are under your PATH.
* APAVplot <br>APAVplot is an R package specifically designed for visualization of PAV analysis. The code and more details please see [here](https://github.com/SJTU-CGM/APAVplot).

## Installation procedures

1. Download the APAV toolbox from github:
```
$ git clone git@github.com:SJTU-CGM/APAV.git
```
 Alternatively, you also could obtain the toolbox in the [APAV](https://cgm.sjtu.edu.cn/APAV/install.html) website and uncompress the APAV toolbox package:
```
$ tar zxvf APAV-v**.tar.gz
```

2. Add `apav` to PATH and add `lib/` to PERL5LIB
```
$ export PATH=$PATH:/path/to/APAV/:
$ export PERL5LIB=$PERL5LIB:/path/to/APAV/lib/:
```

3. Install R package 'APAVplot'
```
## Install "ComplexHeatmap" from BiocManager 
$ conda install r-BiocManager  ## OR: Rscript -e "install.packages('BiocManager')"
$ Rscript -e "BiocManager::install('ComplexHeatmap')" 
## OR install it from bioconda 
$ conda install bioconda::bioconductor-complexheatmap

## Install "APAVplot"
$ conda install r-devtools ## OR: Rscript -e "install.packages('devtools')" 
$ Rscript -e "devtools::install_github('xiaonui/APAVplot')"
```
_You can skip this step if you are going to use APAVplot later to draw plots in the R environment_

4. Test if APAV toolkit is installed successfully
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
        pavSim                  Simulate the size of pan-genome and core-genome from the PAV table
        pavPlotSim              Draw growth curve of genome simulation

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

The main steps can be automatically executed with the `geneBatch`/`generalBatch` command. `geneBatch` command is used for gene region and takes GFF file as input. `generalBatch` command is used for any target region and takes BED file as input.

```
cd ${APAV_PATH}/demo/

## demo1: Some genes on human chromosome 19
apav geneBatch --gff demo1_gene.gff3 --bamdir bam --pheno demo_sample.pheno --fa demo.fa.gz --fam demo1_gene.fam --up_n 10 --down_n 10 --chrl demo1.chrl --out demo1

## demo2: Some proteins on human chromosome 19
apav generalBatch --bed demo2_general.bed --bamdir bam --pheno demo_sample.pheno --fa demo.fa.gz --out demo2

## demo3: Some repeat sequences on human chromosome 19
apav generalBatch --bed demo3_general.bed --bamdir bam --pheno demo_sample.pheno --fa demo.fa.gz --rmele --out demo3
```

The step-by-step procedure is available on the [APAV](https://cgm.sjtu.edu.cn/APAV/start.html) website. The PAV reports for demos can also be viewed on [here](https://cgm.sjtu.edu.cn/APAV/demo.html).

