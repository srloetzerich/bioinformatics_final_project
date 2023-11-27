Run step 1 from within the each of the make_lastz_chains output folders for each species.  
Aggregate these into one directory in step 2,
then run ROAST for step 3

# 0. install multiz:

module load CMake/3.26.3

git clone https://github.com/multiz/multiz.git

cd multiz

make

make install (couldn't get this to run, we shall see what happens)

# 1. NET to MAF – run separately from within each make_last_chains output directory. The sys.argv[1] position = reference species (han in this example)

export PATH=/project/daane/sara/multiz/GenomeAlignmentTools/bin:$PATH

export PATH=/project/daane/shared/make_lastz_chains/GenomeAlignmentTools/kent/bin:$PATH

export PATH=/project/daane/shared/make_lastz_chains/GenomeAlignmentTools/src:$PATH

export PATH=/project/daane/shared/multiz/bin:$PATH

python /home/srloetze/python_script/sbatch_sabine.py --jobname nets2maf --time 25:00 --batchname net --mem 50G --email email@email.com --script "python /home/srloetze/python_script/nets_to_maf.py dr"

#sbatch_sabine.py, original by Jacob Daane

```
#goal to print sbatch files for O2 job submission
import sys
import subprocess

usage = "usage: %prog [options]"
from optparse import OptionParser
parser = OptionParser(usage=usage)

parser.add_option("--threads", help='-t 1', default='1')
parser.add_option("--time", help='-t 12:00', default='12:00')
parser.add_option("--jobname", help='-job-name=sift',default='default')
parser.add_option("--mem", help='--mem=4G', default='4G')
parser.add_option("--email", help='--mail-user=email@gmail.com', default='email@gmail.com')
parser.add_option("--script", help='/n/groups/harris/Sara/sift4g/bin/sift4g  -q all.fa -d /n/groups/harris/Jake/protein_database/uniref90.fasta')
parser.add_option("--batchname", help='tmpbatch',default='tmpbatch')
parser.add_option("--input", help='filein.txt',default=False)
parser.add_option("--output", help='fileout.txt',default=False)

options,args = parser.parse_args()

CPUs = options.threads
T = options.time
name = options.jobname
memory = options.mem
mail = options.email
command = options.script
bname = options.batchname
IN = options.input
OUT = options.output

import os
import os.path
if os.path.isfile(bname):
    rm_cmd = 'rm ' + str(bname)
    subprocess.call(rm_cmd, shell=True)

sys.stdout = open(bname, 'a')

print ('#!/bin/bash')
print ('#SBATCH -J ' + str(name))
print ('#SBATCH -t ' + str(T))
print ('#SBATCH -N 1 -n ' + str(CPUs))
print ('#SBATCH --mem=' + str(memory))
print ('#SBATCH --mail-user=' + str(mail))
print ('#SBATCH --mail-type=ALL')
if IN:
    print (')#SBATCH --input '  + str(IN))
if OUT:
    print ('#SBATCH --output ' + str(OUT))
else:
    print ('#SBATCH -o ' + str(name) + '.0%j')
print ('#SBATCH -e ' + str(name) + '_%j.error')
print ('')
print (command)

sys.stdout.close()

batch_cmd = 'sbatch ' + str(bname)
subprocess.call(batch_cmd,shell=True)

```
```
#updated script using argparse
#requires argparse module in python environment

import sys
import subprocess
import argparse
import os

parser = argparse.ArgumentParser(description="Generate sbatch files for O2 job submission")

parser.add_argument("--threads", help='-t 1', default='1')
parser.add_argument("--time", help='-t 12:00', default='12:00')
parser.add_argument("--jobname", help='-job-name=sift', default='default')
parser.add_argument("--mem", help='--mem=4G', default='4G')
parser.add_argument("--email", help='--mail-user=email@gmail.com', default='email@gmail.com')
parser.add_argument("--script", help='/n/groups/harris/sara/sift4g/bin/sift4g  -q all.fa -d /n/groups/harris/sara/protein_database/uniref90.fasta')
parser.add_argument("--batchname", help='tmpbatch', default='tmpbatch')
parser.add_argument("--input", help='filein.txt', default=False)
parser.add_argument("--output", help='fileout.txt', default=False)

args = parser.parse_args()

CPUs = args.threads
T = args.time
name = args.jobname
memory = args.mem
mail = args.email
command = args.script
bname = args.batchname
IN = args.input
OUT = args.output

if os.path.isfile(bname):
    rm_cmd = 'rm ' + str(bname)
    subprocess.call(rm_cmd, shell=True)

with open(bname, 'a') as sbatch_file:
    sbatch_file.write('#!/bin/bash\n')
    sbatch_file.write('#SBATCH -J ' + str(name) + '\n')
    sbatch_file.write('#SBATCH -t ' + str(T) + '\n')
    sbatch_file.write('#SBATCH -N 1 -n ' + str(CPUs) + '\n')
    sbatch_file.write('#SBATCH --mem=' + str(memory) + '\n')
    sbatch_file.write('#SBATCH --mail-user=' + str(mail) + '\n')
    sbatch_file.write('#SBATCH --mail-type=ALL\n')
    if IN:
        sbatch_file.write('#SBATCH --input ' + str(IN) + '\n')
    if OUT:
        sbatch_file.write('#SBATCH --output ' + str(OUT) + '\n')
    else:
        sbatch_file.write('#SBATCH -o ' + str(name) + '.0%j\n')
    sbatch_file.write('#SBATCH -e ' + str(name) + '_%j.error\n\n')
    sbatch_file.write(command)

batch_cmd = 'sbatch ' + str(bname)
subprocess.call(batch_cmd, shell=True)

```
#nets_to_maf.py
 
 ```
import sys
import subprocess
import glob

# Identify species acronyms
reference = sys.argv[1] # acronym used throughout the files:  han, cgo, etc.
twobit_files = glob.glob('*2bit')
query = ''
for infile in twobit_files:
    spec = infile.split('.')[0]
    if spec != reference:
        query = spec

# 1. chainNet (needs memory)
chainNet_cmd = 'chainNet ' + str(reference) + '.' + str(query) + '.allfilled.chain.gz ' + str(reference) + '.chrom.sizes ' + str(query) + '.chrom.sizes ' + str(reference) + '.net ' + str(quer$
subprocess.call(chainNet_cmd, shell=True)

# 2. netSyntenic (fast)
netSyntenic_cmd = 'netSyntenic ' + str(reference) + '.net ' + str(reference) + '_ns.net'
subprocess.call(netSyntenic_cmd, shell=True)

# 3. netFilter (fast)
netFilter_cmd = 'NetFilterNonNested.perl -doScoreFilter -keepSynNetsWithScore 5000 -keepInvNetsWithScore 5000 -minScore1 100000  ' + str(reference) + '_ns.net > ' + str(reference) + '.filtere$
subprocess.call(netFilter_cmd, shell=True)

# 4. NET to AXT -- 4 min but 50Gb needed
nettoAxt_cmd = 'netToAxt ' + str(reference) + '.filtered.net ' + str(reference) + '.' + str(query) + '.allfilled.chain.gz ' + str(reference) + '.2bit ' + str(query) + '.2bit stdout | axtSort $
subprocess.call(nettoAxt_cmd, shell=True)

# 5. AXT to MAF -- fast, low mem
axtToMaf_cmd = 'axtToMaf ' + str(reference) + '-' + str(query) + '-filtered.axt ' + str(reference) + '.chrom.sizes ' + str(query) + '.chrom.sizes ' + str(reference) + '-' + str(query) + '-fil$
subprocess.call(axtToMaf_cmd, shell=True)

# 6.  single cov
single_cov_cmd = 'single_cov2 ' + str(reference) + '-' + str(query) + '-filtered.maf > ' + str(reference) + '.' + str(query) + '.sing.maf'
subprocess.call(single_cov_cmd, shell=True)
 
 ```

# 2. Aggregate sing.maf files output from #1 into one directory for all species for multiz. Rename as ref.query.sing.maf (e.g. : han.dlo.sing.maf)

` done manually.... unless file names are consistent across directories `

#move all sing.maf files to one directory

```
#!/bin/bash

# Define the source and destination directories
source_dir="/path/to/source"
destination_dir="/path/to/destination/sing_maf"

# Create the destination directory if it doesn't exist
mkdir -p "$destination_dir"

# Loop through directories ending in "chain"
for dir in "${source_dir}"/*chain; do
    if [ -d "$dir" ]; then
        # Loop through files ending in "sing.maf"
        for file in "$dir"/*sing.maf; do
            if [ -f "$file" ]; then
                # Copy the file to the destination directory
                cp "$file" "${destination_dir}/$(basename "$file")"
            fi
        done
    fi
done


```
# 3. Run multiz on all files to aggregate into one multiple sequence alignment 

(run from within directory containing all the sing.maf files from step 2). 
T= tmp directory that you want multiz to use, I recommend a tmp directory within the working directory. 
E= reference species. 
Newick tree is required formatted as shown (no commas). 
*.sing.maf tells multiz to look for files ending in .sing.maf. 
nonicefishes.maf is the output maf

python /home/srloetze/python_script/sbatch_sabine.py --jobname multiz --time 36:00:00 --batchname net --mem 50G --email email@email.com --script "roast + X=0 T=/project/daane/sara/multiz/sing_maf/tmp E=dr '(((((((((((li cv) gg) ci) mam) dr) (dt1 dd)) om) pc) ma) (my xyt)))' *.sing.mafq cypriniforms.maf"

```
#ran in sing_maf/ to get it to work, not multiz/

#tree
(((((((((((li,cv),gg),ci),mam),dr),(dt1,dd)),om),pc),ma),(my,xyt)));
 
 ```
