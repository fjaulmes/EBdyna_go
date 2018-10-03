echo "********************************************"
echo "First make sure that the necessary mexa64 files are in the main directory"
cp ./ba_interp2/*.mexa64 ./
cp ./ba_interp3/*.mexa64 ./
cp ./sharedmatrix/*.mexa64 ./
ls *.mexa64

echo "*********** cleaning old log files *********"
rm ~/java.log.*
rm ./$1*.o*
rm ./$1*.e*

echo "********************************************"
echo "Specified simulation ID: $1"
echo "Starting compilation of EBdyna_go (load_RAM and G_eq), please wait..."


echo "********************************************"
echo "removing any previous files in specified output"
rm ./output_$1/*.mat 

echo "Starting compilation load_RAM, please wait..."


cp -r ../01\ general_functions/* ./
module load matlab/2016a
# mcc -m load_RAM.m
rm ./load_RAM
file="./load_RAM"
if [ -f "$file" ]
then
	echo "Using previously compiled load_RAM."
else
	echo "$file not compiled yet."
	mcc -m load_RAM.m
	echo "Succesfully compiled load_RAM."
fi


echo "Starting compilation G_eq, please wait..."

mcc -R -singleCompThread -m G_eq.m
rm -f run_G_eq.sh
rm -f run_load_RAM.sh
rm -f *.txt
rm -f *.log
echo "Succesfully compiled G_eq"

echo "Coverting compiled file to ID: $1"
mkdir -p ./output_$1
cp ./G_eq ./G_eq_$1
cp ./load_RAM ./load_RAM_$1

echo "Resetting the pool for ID: $1"
./reset_pool_i.sh $1

echo "EBdyna_go compiled succesfully!"
echo "********************************************"

# rm ~/java.log.*
