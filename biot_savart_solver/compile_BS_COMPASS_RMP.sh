mkdir -p ./output_RMP
cp -r  ../general_functions/*.m ./
module load matlab/2016a
mcc -R -singleCompThread -m BS_COMPASS_RMP_main.m
rm -f run_BS_COMPASS_RMP_main.sh
rm -f *.txt
rm -f *.log

./reset_pool_RMP.sh