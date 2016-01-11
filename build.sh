
make -C src/ clean
make -C src/

EXPERIMENTS_DIR="data"
OUT_DIR_NAME=`hostname`_`date +%F`_`date +%R | xargs`

OUTPUT_DIRECTORY="${EXPERIMENTS_DIR}/${OUT_DIR_NAME}"
mkdir -p $OUTPUT_DIRECTORY
EXPR_FILE="${OUTPUT_DIRECTORY}/experiments.csv"
RESULT_FILE="${OUTPUT_DIRECTORY}/experiments_results.csv"
FINAL_RESULTS_DIR="${EXPERIMENTS_DIR}/last"

NB_THREADS=8

# generate the DoE
Rscript ./scripts/gen_experiments.R $EXPR_FILE

# run the experiments
python ./scripts/run_experiments.py $EXPR_FILE $RESULT_FILE

# analyze the results
Rscript ./scripts/analysis.R $RESULT_FILE $OUTPUT_DIRECTORY $NB_THREADS

# move the generated pdf to the output directory
mv Rplots.pdf $OUTPUT_DIRECTORY

# create a symbolic link to the last experiment results
rm -f $FINAL_RESULTS_DIR
ln -s $OUT_DIR_NAME $FINAL_RESULTS_DIR
