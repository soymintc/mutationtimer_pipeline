[ $# -ne 4 ] && { echo -e "\nUsage: $0 <base.dir> <sample_id> <maf> <remixtpp>\n" 1>&2; exit 1; }
indir=$1
sample_id=$2
maf=$3
remixtpp=$4
[ ! -f $maf ] && { echo "ERROR: $maf does not exist"; exit 1; }

base_dir=$(realpath $indir)
[ ! -d $base_dir ] && { echo "ERROR: $base_dir does not exist"; exit 1; }

CLUSTER_CMD=("bsub -n {threads} -R {cluster.resources} -M {cluster.memory} -o {cluster.output} -J {cluster.name} -W {cluster.time}")

results_dir=$base_dir/results
log_dir=$base_dir/log
cluster_yaml=cluster.yaml

[ ! -d $results_dir ] && { echo "LOG: $results_dir does not exist" 1>&2; mkdir -p $results_dir; }
[ ! -d $log_dir ] && { echo "LOG: $log_dir does not exist" 1>&2; mkdir -p $log_dir; }
[ ! -f $cluster_yaml ] && { echo "ERROR: $cluster_yaml does not exist" 1>&2; exit 1; }

cmd="snakemake --config"
cmd="$cmd results_dir=$results_dir"
cmd="$cmd log_dir=$log_dir"
cmd="$cmd sample_id=$sample_id"
cmd="$cmd maf=$maf"
cmd="$cmd remixtpp=$remixtpp"
cmd="$cmd --jobs 10"
cmd="$cmd --use-singularity"
cmd="$cmd --singularity-args \"--bind /juno\""
cmd="$cmd --skip-script-cleanup"
cmd="$cmd --cluster-config $cluster_yaml"
cmd="$cmd --cluster \"${CLUSTER_CMD}\""
cmd="$cmd --cluster-cancel bkill"
cmd="$cmd --snakefile mutationtimer.smk"
cmd="$cmd --rerun-incomplete" #--dry-run" #--dag

echo $cmd
eval $cmd
