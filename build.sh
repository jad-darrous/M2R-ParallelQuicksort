
# global variables
EXPERIMENTS_DIR="data"
OUT_DIR_NAME=`hostname`_`date +%F`_`date +%R | xargs`

OUTPUT_DIRECTORY="${EXPERIMENTS_DIR}/${OUT_DIR_NAME}"
mkdir -p $OUTPUT_DIRECTORY

EXPR_FILE_1="${OUTPUT_DIRECTORY}/experiments_small.csv"
EXPR_FILE_2="${OUTPUT_DIRECTORY}/experiments_medium.csv"
EXPR_FILE_3="${OUTPUT_DIRECTORY}/experiments_big.csv"

RESULT_FILE_1="${OUTPUT_DIRECTORY}/experiments_small_results.csv"
RESULT_FILE_2="${OUTPUT_DIRECTORY}/experiments_medium_results.csv"
RESULT_FILE_3="${OUTPUT_DIRECTORY}/experiments_big_results.csv"

FINAL_RESULTS_DIR="${EXPERIMENTS_DIR}/last"

NB_THREADS=8


# generate the DoE
Rscript ./scripts/gen_experiments.R $EXPR_FILE_1 $EXPR_FILE_2 $EXPR_FILE_3


# make sure that the binaries are up-to-date
make -C src/ clean
make -C src/ ARGS="-DALL_IMPL"

# run the experiments
python ./scripts/run_experiments.py $EXPR_FILE_1 $RESULT_FILE_1

# run the experiments
python ./scripts/run_experiments.py $EXPR_FILE_2 $RESULT_FILE_2

# make sure that the binaries are up-to-date
make -C src/ clean
make -C src/ ARGS=""

# run the second experiments
python ./scripts/run_experiments.py $EXPR_FILE_3 $RESULT_FILE_3


# generate the report
if [ -x "$(command -v pandoc)" ]; then
  # requires pandoc
  Rscript -e "library(rmarkdown); render('report.Rmd')"
  mv report.html $OUTPUT_DIRECTORY
else
  echo "pandoc is required to generate the report."
  echo "generating a simple version instead"
  # generate md file
  Rscript -e "library(knitr); knit('report.Rmd')"
  mv report.md figure $OUTPUT_DIRECTORY
fi

# create a symbolic link to the last experiment results
rm -f $FINAL_RESULTS_DIR
ln -s $OUT_DIR_NAME $FINAL_RESULTS_DIR
