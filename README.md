# TMDB_grid

Code for TMDB data analysis using CERN grid. The code generates a lot of histograms for TMDB performance analysis

Instructions:
1. Check if GCC and ROOT are installed
2. Execute: make clean; make
3. Running locally: ./TMDB $(find your_dataset_directory -name "*.root" -printf '%p,') output_file
4. Running on grid: edit submit.py with your options and run
