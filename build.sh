
make -C src/ clean
make -C src/

OUTPUT_DIRECTORY=data/`hostname`_`date +%F`_`date +%R | xargs`
mkdir -p $OUTPUT_DIRECTORY
EXPR_FILE="${OUTPUT_DIRECTORY}/experiments.csv"
RESULT_FILE="${OUTPUT_DIRECTORY}/experiments_results.csv"
NB_THREADS=8

Rscript ./scripts/gen_expr.R $EXPR_FILE
python ./scripts/run_experiments.py $EXPR_FILE $RESULT_FILE
Rscript ./scripts/analysis.R $RESULT_FILE $OUTPUT_DIRECTORY $NB_THREADS
mv Rplots.pdf $OUTPUT_DIRECTORY
