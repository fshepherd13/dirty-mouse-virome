## Dissemination analysis workflow

#### The workflow for identifying variants that disseminate from the small intestine to the liver works slightly differently from the workflow described in ~/transmission_bottlenecks/README.md. This is because instead of calling variants against a single pet store consensus, each mouse needs its liver reads compared to its own small intestines. Therefore, consensus sequences need to be generated for all mice.

### Generating consensus sequences for small intestine reads.
#### Use the snakemake file under ./si_consensus_generation/ to do this step. Similarly to the transmission bottleneck snakefile, this creates a consensus by aligning reads generated from small intestine virus amplicons. 

#### First, create a directory to hold the raw reads from both the small intestine and liver amplicon seq.
```
mkdir raw_reads
```

#### Move the raw fastq.gz reads to this folder. Unlike the transmission bottlenecks, it isn't necessary to group reads by cage for analysis. But your raw reads do need to have a unique name per animal so the snakefile can find the correct SI consensus sequence to call liver variants. For example a pair of corresponding reads should be named like this:
#### `M1_AK04-SI-astro-rdrp-D2_A_R1.fastq.gz` <- SI reads from mouse 1, replicate A
#### `M1_AK04-SI-astro-rdrp-D2_A_R2.fastq.gz` <- " " "

#### `M1_AK04-LIV-astro-rdrp-D2_A_R1.fastq.gz` <-  Liver reads from mouse 1 replicate A
#### `M1_AK04-LIV-astro-rdrp-D2_A_R2.fastq.gz` <-  " " "

#### Notice that everything in the name is the same, except for the "LIV" and "SI" designations. 

#### Navigate to the directory `si_consensus_generation`
```
$ cd si_consensus_generation/
```
#### Run the snakemake file in this directory using a command similar to the following, which will generate a consensus sequence for each of the small intestine reads, changing the sample ID's to match your experiment:
```
$ snakemake --cores 3 out/consensus/M{1,3,4,5,10,12}_AK04-SI-astro-rdrp-D2.fa
```

#### You can also create FASTQC files for raw and trimmed reads with this code:
```
$ snakemake --cores 2 quality/fastqc_trimmed/M{1,3,4,5,10,12}_AK04-SI-astro-rdrp-D2_{A,B}_R{1,2}_trimmed_fastqc.html
$ snakemake --cores 2 quality/fastqc_raw/M{1,3,4,5,10,12}_AK04-SI-astro-rdrp-D2_{A,B}_R{1,2}_fastqc.html
```

#### Inspect the assembled consensus sequences to make sure they are the correct size and bind with the PCR primers. I do this in Geneious. I also annotate the correct open reading frame of the amplicon and export a .gff3 annotation file from Geneious so that the variant call pipeline can identify amino acid changes. I always export the gff file to a new folder (that I create with mkdir command) ./out/gff.

### Call variants in liver and small intestine reads against the small intestine consensus sequence.
#### In order to perform downstream analyses, you not only need to call variants that occur in the liver sequences, you also need to call variants within the SI sequences. This allows for graphing the frequency of variants that occur during dissemination. 

#### Navigate to the `variant_calling` directory.
```
$ cd ../variant_calling/
```

#### Run the snakefile in this directory. I first two commands from the command line, one to call variants in the LIV reads:
```
$ snakemake --cores 3 variant_calling/variants/M{1,3,4,5,10,12}_AK04-LIV-astro-rdrp-D2.filtered.tsv
```
#### ...and one to call variants in the SI reads:
```
$ snakemake --cores 3 variant_calling/variants/M{1,3,4,5,10,12}_AK04-SI-astro-rdrp-D2.filtered.tsv
```

#### The combined variant call files located in ./out/variants/ need to be aggregated, which I do with a python script. The script requires an outside file of metadata (see metadata.csv file for example) to add experimental info onto that gets tacked on to the aggregated csv file. The python script also calculates an average variant frequency for the two pcr replicates. An example code to run follows (see python script as well for general usage):
```
python3 /combine_variant_files.py example_metadata.csv ./out/variants/ combined_filtered_dissem_variants.csv
```

#### Variants can now be graphed using code in `~/figures/`.

