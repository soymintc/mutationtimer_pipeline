[ $# -ne 2 ] && { echo -e "\nUsage: $0 <base.dir> <sample_id>\n" 1>&2; exit 1; }
indir=$1
sample_id=$2

base_dir=$(realpath $indir)
[ ! -d $base_dir ] && { echo "LOG: $base_dir does not exist"; mkdir -p $base_dir; }

CLUSTER_CMD=("bsub -n {threads} -R {cluster.resources} -M {cluster.memory} -o {cluster.output} -J {cluster.name} -W {cluster.time}")

results_dir=$base_dir/results
intermediate_dir=$base_dir/tmp
log_dir=$base_dir/log
cluster_yaml=config/cluster.yaml
pipeline_config=config/pipeline.yaml

[ ! -d $results_dir ] && { echo "LOG: $results_dir does not exist" 1>&2; mkdir -p $results_dir; }
[ ! -d $intermediate_dir ] && { echo "LOG: $intermediate_dir does not exist" 1>&2; mkdir -p $intermediate_dir; }
[ ! -d $log_dir ] && { echo "LOG: $log_dir does not exist" 1>&2; mkdir -p $log_dir; }
[ ! -f $cluster_yaml ] && { echo "ERROR: $cluster_yaml does not exist" 1>&2; exit 1; }

cmd="snakemake --config"
cmd="$cmd results_dir=$results_dir"
cmd="$cmd intermediate_dir=$intermediate_dir"
cmd="$cmd log_dir=$log_dir"
cmd="$cmd sample_id=$sample_id"
cmd="$cmd --configfile $pipeline_config"
cmd="$cmd --jobs 10"
cmd="$cmd --use-singularity"
cmd="$cmd --singularity-args \"--bind /juno\""
cmd="$cmd --skip-script-cleanup"
cmd="$cmd --cluster-config $cluster_yaml"
cmd="$cmd --cluster \"${CLUSTER_CMD}\""
cmd="$cmd --cluster-cancel bkill"
cmd="$cmd --rerun-incomplete" #--dry-run" #--dag

echo $cmd
eval $cmd
