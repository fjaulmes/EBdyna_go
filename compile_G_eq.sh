echo "********************************************"
export EBdyna_rundir=d3d_184342_2500_15R_01
echo "Coying 3 main files from EBdyna_rundir=" $EBdyna_rundir
cp ./$EBdyna_rundir/G_eq.m ./
cp ./$EBdyna_rundir/define_basic_parameters.m ./
cp ./$EBdyna_rundir/sim_parameters.m ./
echo "First make sure that the necessary mexa64 files are in the main directory"
cp -r ./EBdyna/G_eq/init_data/*.m ./
cp ./EBdyna/G_eq/general_functions/* ./
cp -r ./EBdyna/G_eq/trajectories/*.m ./
cp ./EBdyna/G_eq/GT_eq.m ./
#cp -r ../data_common/sigma_cx*.mat ../data_tokamak
#cp -r ../data_common/physics_constants.mat ../data_tokamak
mkdir ./$EBdyna_rundir/output_singlefile


echo "*********** cleaning old log files *********"
rm ~/java.log.*
rm ./*.out
rm ./*.e*
rm ./G_eq

wait
echo "********************************************"
echo "Specified simulation ID: $1"
echo "Starting compilation of EBdyna_go (G_eq), please wait..."


echo "********************************************"

# echo "Starting compilation load_RAM, please wait..."



module load MATLAB/R2021a
# mcc -m load_RAM.m



echo "Starting compilation G_eq, please wait..."

mcc -R -singleCompThread -m G_eq.m
rm -f run_G_eq.sh
rm -f run_load_RAM.sh
rm -f *.txt
rm -f *.log
rm -f *.m
rm -f *.mexa64
echo "Succesfully compiled G_eq"

echo "EBdyna_go compiled succesfully!"
echo "********************************************"

# rm ~/java.log.*

