# **Convergent Gene Loss & Duplication Analysis of Cypriniforms**

For this, I'll be starting with 'TOGA'-- [Tool to infer Orthologs from Genome Alignments.](https://github.com/hillerlab/TOGA#readme) TOGA has an established pipeline consisting of:
*Repeat Modeler
*Repeat Masker
*Make LastZ Chains
*Run TOGA #TOGA Party

My lab and I have also come up with some additionaly analysis for our TOGA outputs which will be shared later. First, I will download genomes for my species of interest.

## Downloading & Scrubbing Genomes
*my genomes were sourced from [NCBI](https://www.ncbi.nlm.nih.gov/) but many can also be found at [Ensembl](https://ensemblgenomes.org/).
*I would download with the command:

```bash
curl -OJX GET "linkaddress" #substituting the link address of the genome of your choice
```
Then I would unzip using 'unzip filename'. 

Next, I would remove all periods and spaces since this is a known hiccup when using the TOGA pipeline. 

```bash
sed 's/ .*//g' input_file.fa > output_file.fa #removes spaces
sed 's/\./_/g' input_file.fa > output_file.fa  #removes periods and replaces them with '_'
```

In my next section, I begin actually running jobs for this pipeline-- almost all of my jobs have the same batch script heading which I'll share here and skip later to avoid redundancy.

### Repeat Modeler
This package allows for the identification and classification of repeated elements in a genomic sequence of interest. 

```bash
#!/bin/bash

#SBATCH -J rmod #jobname

#SBATCH -o rmod_species.%j #jobname

#SBATCH -t 335:00:00 #how long

#SBATCH -N 1 -n  #how many CPUs

#SBATCH --mem=10G #how much memory

#SBATCH --mail-user=email #where to send job email with begin/end

#SBATCH --mail-type=ALL #clarify what is being mailed

module load RepeatModeler/2.0.4-foss-2022a #only works if module is already installed, install instructions on the TOGA github

BuildDatabase -name genomenamehere path/to/genome

RepeatModeler -database dt2 -threads 28 -LTRStruct
```
to run: 'sbatch /path/to/batchscript'

*Of note, genomes with poor quality will struggle with this modeling. I found it worked well to model the trouble genomes to the reference genome instead where necessary.
#### Repeat Masker
*This step looks at the repeated portions of the genomic sequence and 'masks' highly repetetive segments with 'N', hard masking or a lowercase version of their original nucleotide 'agct', soft masking. Both are viable options, my lab group found that soft masking gave better analysis overall.

To run:

```bash
module load RepeatMasker/4.1.5nano-foss-2022a

RepeatMasker -lib repeatmodeloutput-familes.fa -e rmblast -xsmall -speciesname 28 -dir /path/to/masked/output -s /path/to/original/fasta.fa
```
##### Make LastZ Chains

There are some essential python modules for this portion:
```bash
module load Python/3.9.6-GCCcore-11.2.0
cd /path/to/location/of/make_lastz_chains
pip3 install -r requirements.txt
```
Before running any of these jobs, be sure to add these to your paht and ensure python **is loaded**.

```bash
PATH=$PATH:/path/to/shared/lastz-distrib/bin;export PATH

PATH=$PATH:/path/to/shared/make_lastz_chains/GenomeAlignmentTools/kent/bin; export PATH

export KENTSRC_DIR=/path/to/shared/make_lastz_chains/GenomeAlignmentTools/kent/src/

module add Nextflow/22.10.6 

python /path/to/make_lastz_chains/make_chains.py referencegenome querygenome  /path/to/maskedreferencegenome.fa path/to/maskedquerygenome.fa --project_dir nameofdirectory --executor slurm --executor_queuesize 200 --seq1_chunk 50000000 --seq2_chunk 10000000 
```
Sometimes it was necessary to change chunk sizes to get larger genomes to run in a tolerable timeframe. Note, this kicks off many jobs and takes a long time. It is easily the most difficult of all the steps combined.
*outputs will include one .2bitfile for each query and reference species as well as an allfilled.chain.gz; these are needed for the next step
*you will also need a bed file  and isoform file for your reference genome, this can be downloaded at either NCBI or ensemble. It is recommended to be consistent, source this from your genome source location.

##### TOGA Party
Before running: 
```bash
module load Python/3.10.6-GCCcore-11.3.0

pip3 install py_nf --user

pip3 install -r requirements.txt --user  (must be in TOGA location)
```
to run:
```bash
module add Nextflow/22.10.6

export PATH=/home/srloetze/.local/lib/python3.7/site-packages:$PATH #use whichever works for you!

export PATH=/home/srloetze/.local/lib/python3.9/site-packages:$PATH #use whichever works for you!

python /path/to/shared/TOGA/toga.py /path/to/refgenome.querygenome.allfilled.chain.gz /path/to/referencegenomebedfile.bed /path/to/referencegenome.2bit /pathtoquerygenome.2bit --pn querygenomename -i /path/to/referencegenomeisoform.txt --nc /path/to/TOGA/nextflow_config_files --cesar_buckets 2,10,20,50
```

This will give you many different output files that can be seen on the Hiller Lab Github for TOGA. I am most interested in the file 'loss_summ_data.tsv'. I want to collect them for all of my species and annotate them with their gene names and description. I also want to see the percentage that is intact/partially intact per this analysis for quality control.

#TOGA Analysis

In order to pull out loss_sum or other file and put into a single output folder:

```bash
#!/bin/bash

source_directory="/path/to/source_dir/" # Set the source directory where the series of directories are located

output_directory="/path/to/output_dir/" # Set the output directory where the copied files will be placed

mkdir -p "$output_directory" # Ensure the output directory exists

for dir in "$source_directory"/*; do # Iterate through the source directories
    if [ -d "$dir" ]; then       
        directory_name=$(basename "$dir") # Extract the last three letters from the directory name
        suffix="${directory_name: -3}"     # Check if the "loss_sum_data.tsv" file exists in the directory
        if [ -e "$dir/loss_sum_data.tsv" ]; then 
            cp "$dir/loss_sum_data.tsv" "$output_directory/$suffix.tsv" # Copy the file to the output directory with the desired name
            echo "Copied $dir/loss_sum_data.tsv to $output_directory/$suffix.tsv"
        else
            echo "No 'loss_sum_data.tsv' found in $dir"
        fi
    fi
done
```

Before running the following script:

module load Python/3.9.6-GCCcore-11.2.0

pip install numpy --user

python /path/to/sbatch_sabine.py --jobname nameofyourjob --time 25:00 --batchname nameofyourjobbatch --mem 50G --email youremailhere@gmail.com --script "python /path/to/test_toga_composition.py < loss_summ_data.tsv"

To evaluate composition of the loss_sum_data.tsv:

*sbatch_sabine.py is a script that writes my batch script for me with this input
*test_toga_composition.py listed below, originally written by the lab PI, Jacob Daane.

```bash
import sys
import csv
import re
from collections import defaultdict
import numpy as np

total = 0
datDict = defaultdict(int)
for f in csv.reader(sys.stdin,csv.excel_tab):
    if f[0] == 'GENE':
        dat = f[2]
        datDict[dat] += 1
        total += 1

for item in datDict:
    percentage = float(datDict[item])/float(total)
    print (str(item) + '\t' + str(datDict[item]) + '\t' + str(percentage))

intact_group = datDict['PI'] + datDict['I']
intact_percentage = float(intact_group)/total
print('I+PI' + '\t' + str(intact_group) + '\t' + str(intact_percentage))

loss_group = datDict['L'] + datDict['UL']
loss_percentage = float(loss_group)/total
print('L+UL' + '\t' + str(loss_group) + '\t' + str(loss_percentage))
```
Now we know what we are working with, we can annotate our gene lists! 

#Annotate with Correct Gene(s)

You will want to go to the Bio Mart in Ensembl. Once there you will drop down "choose database" and select Ensembl genes. This will prompt you to choose a species that you want to retrieve information from.

After you get to this step you will want to select attributes, and drop down the gene option. There is plenty of options to choose from and you will click on the attributes you are interested in. From there you will click on results and the file will download.

Then combine all information using the following python script

#Annotation Script

```bash
#!/bin/bash

output_file="output.tsv"
annotation_file="mart_export.txt"
output_with_descriptions="output_with_descriptions.tsv"

# Create a temporary file to store the merged data
temp_file=$(mktemp)

# Use awk to join the files based on the gene ID
awk -F'\t' 'NR == FNR { gene_desc[$1] = $0; next } FNR == 1 { print $0 "\tGene_Description"; next } { if ($1 in gene_desc) { print $0 "\t" gene_desc[$1] } else { print $0 "\tNA" } }' "$annotation_file" "$output_file" > "$temp_file"

# Rename the temporary file to the desired output file
mv "$temp_file" "$output_with_descriptions"

echo "Gene descriptions added and saved as $output_with_descriptions"
```     
