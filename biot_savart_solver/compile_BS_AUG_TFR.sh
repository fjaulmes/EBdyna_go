mkdir -p ./output_TF
cp -r  ../01\ general_functions/*.m ./
module load matlab/2016a
mcc -R -singleCompThread -m BS_AUG_TFR_main.m
rm -f run_BS_AUG_TFR_main.sh
rm -f *.txt
rm -f *.log

./reset_pool_TFR.sh