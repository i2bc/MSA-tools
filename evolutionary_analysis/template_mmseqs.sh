#! /bin/sh

#PBS -N colabsearch
#PBS -q alphafold
#PBS -o $LOGDIR/colabsearch.output.log
#PBS -e $LOGDIR/colabsearch.error.log
#PBS -l mem=8gb
#PBS -l ncpus=1
#PBS -l walltime=02:00:00
#PBS -l node_group=ALIGN

# VARIABLES:
# LOGDIR => where pbs log files will be saved
# OUTDIR => where the msa will be generated
# PROGRESSLOG => name and path to the job progress log file (for printing on the webserver)
# INPUT => name and path to the fasta input file
# JOBSCRIPT => path to jobcmd.sh for r4s pipeline
# JOBIDFILE => where cluster jobid is saved

set -x
echo Host $(hostname)
echo Jobid $PBS_JOBID

#======= USER COMMANDS =======

echo Start job: $(date)
mkdir -p $OUTDIR
cd $OUTDIR

# write log to file
echo "[$(date '+%F %T')] Loading database..." >> $PROGRESSLOG
DB=/database/db/uniref30_2202/uniref30_2202_db.idx

# check for previous DBs
ISVM=$(ps aux | grep vm | grep -v grep | wc -l)

if [[ $ISVM -eq 1 ]]; then
    # a DB is loaded, check if it's the right one
    ISDB=$(ps aux | grep $DB | grep -v grep | wc -l)
    if [[ $ISDB -eq 0 ]]; then
        # the wrong DB is in cache, it needs to be unloaded first
        OLDDB=$(ps aux | grep vm | grep -v grep | awk '{print $NF}')
        sudo /usr/local/bin/mmseqs_unload -f $OLDDB
    fi
fi

# load the DB
sudo /usr/local/bin/mmseqs_load -f $DB

# write log to file
echo "[$(date '+%F %T')] Running step 1 - colabfold_search msa generation..." >> $PROGRESSLOG

# run colabfold_search
singularity exec --containall --env TMPDIR=$TMPDIR -B $PWD,/store,/data,/opt/alphafold/scripts:/scripts,/opt/alphafold/parameters:/parameters,/database /opt/alphafold/colabfold_container_v15.sif colabfold_search     $INPUT uniref30_2202_db $OUTDIR --db1 /database/db/uniref30_2202/uniref30_2202_db --max-accept 1000000     --use-env 0 --use-templates 0 --threads 1     --db-load-mode 2 -s 8 --filter 1 --expand-eval inf     --align-eval 10 --diff 3000 --qsc 20.0     --mmseqs /opt/MMseqs2/build/bin/mmseqs --db2 . --db3 colabfold_envdb_202108_db

# check for colabfold_search output
myarray=(`find $OUTDIR -maxdepth 1 -name "*.a3m"`)

if [ ${#myarray[@]} -eq 0 ]; then 
    echo "Error: colabfold_search (MMseqs2) step ended with an error (mmseqs outputs weren't generated in $(basename $OUTDIR))" >> $PROGRESSLOG
    exit 1
fi

# write log to file
echo "[$(date '+%F %T')] ...output ok" >> $PROGRESSLOG 

# rename the output a3m after colabfold_search
cd $OUTDIR
echo "[$(date '+%F %T')] ...renaming MSA output files" >> $PROGRESSLOG

for file in $(ls *.a3m); do
    filename=$(python3 <<EOF

import re
input = open("${file}", 'r')
while input:
    line=input.readline()
    if line.startswith('>'):
        name=re.split("\s+",line.replace('>','').strip())[0]
        print(name)
        break
EOF
)

    # if output name already exists, we add a number at the end

    i=1
    while [[ -e ${filename}".a3m" ]]; do
        i=$((i+1))
    done

    if [[ $i -eq 1 ]]; then
        mv $file ${filename}".a3m"
    else
        mv $file ${filename}"_"${i}".a3m"
    fi
done

# write log to file
echo "[$(date '+%F %T')] ...done!" >> $PROGRESSLOG 

if [[ $JOBSCRIPT ]]; then
  qsub -l mem=8Gb -N msa_tools_srv -q common -j oe -o $LOGDIR/ $JOBSCRIPT > $JOBIDFILE
fi
