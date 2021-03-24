
## Workflow description for calling variants in amplicons against a pet store mouse (reservoir host) consensus sequence
#### For this workflow, I've made two folders, each with its own snakemake file. One serves the purpose of creating consensus sequences by amplicon seq reads to a reference sequence. The other compares reads from all mice within a cage to that consensus sequence to call variants. I've been placing raw reads from all mice within a cage within the "consensus_calling" folder.

### Generating a consensus sequence from a pet store mouse

#### First, navigate to the directory called "consensus_sequences"
```
cd consensus_calling
```

#### Create a directory to hold the raw amplicon sequencing reads
```
mkdir raw_reads
```
#### Move the raw reads that belong to one cage of mice to this directory.

#### The snakefile in this directory contains code that trims adapters and low quality reads from raw fastq files, generates FASTQC reports on both raw and trimmed reads (for reference) and generates a consensus sequence by aligning trimmed reads to a specified reference sequence (either Genbank reference sequences, or in the case of divergent sequences, a denovo assembled contig from the raw reads).

#### All mice reads need reads trimmed, but only the pet store mice needs a consensus sequence built. So I run the following commands from the terminal to perform this:
```
$ snakemake --cores 2 trimmed_reads/M{1,3,4}-SI-astro-capsid-E2_{A,B}_R{1,2}_trimmed.fastq #Trims raw read files
$ snakemake --cores 2 quality/fastqc_trimmed/M{1,3,4}-SI-astro-capsid-E2_{A,B}_R{1,2}_trimmed_fastqc.{html,zip} #Runs FASTQC on trimmed reads
$ snakemake --cores 2 quality/fastqc_raw/M{1,3,4}-SI-astro-capsid-E2_{A,B}_R{1,2}_fastqc.{html,zip} #Runs FASTQC on raw reads
$ snakemake --cores 1 out/consensus/M4-SI-astro-capsid-E2.fa #Creates consensus sequence for just the pet store mouse reads
```

#### It is *essential* to check that the consensus sequence assembled here matches the expected amplicon size, doesn't have any weird gaps, and that the primer sequences bind to the ends. I export the sequence to Geneious to check this. I also annotate the correct open reading frame of the amplicon and export a .gff3 annotation file from Geneious so that the variant call pipeline can identify amino acid changes. I always export the gff file to the ./out/index file that the above snakefile will generate for the consensus sequence. 

### Calling variants in trimmed reads against the pet store mouse consensus

#### First, navigate to the variant_call folder to use the correct snakefile.
```
cd ../variant_callings
```

#### The inputs for this snakefile are the trimmed reads generated in above that are stored in the consensus_calling folder. It also requires input for where the assembled petstore consensus sequence is (fasta, index, and gff files). See the comments within the snakefile for more details.

#### There is no need to specify individual target rules or sample names for this, as variants need to be called using every animal's reads. 
```
$ snakemake --cores 2
```

#### The combined variant call files located in ./out/variants/ need to be aggregated, which I do with a python script. The script requires an outside file of metadata (see metadata.csv file for example) to add experimental info onto that gets tacked on to the aggregated csv file. The python script also calculates an average variant frequency for the two pcr replicates. An example code to run follows (see python script as well for general usage):
```
python3 ./combine_variant_files.py metadata.csv ./out/variants/ combined_filtered_variants.tsv
```