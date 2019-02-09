python3 xrun.py 
python process_data.py 
mpirun -n 10 python3 Postprocessing/optimize_driver.py
python process_data.py optim
rm summary.dat 
rm optim_dvars_summary.dat 
rm optim_summary.dat 
