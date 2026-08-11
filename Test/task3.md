### Do simpleISP
1. first you have to make json parameter files for each sub-directories.
   referencing ~/Workspace/simpleISP/example_params.json as sample form.
   referencing parameters.txt 
   infering width, height, bits_per_pixel of params.json from parameters.txt
   raw_path: each raw file path from sub-directories
   output_path: getting each raw file name but change suffix as 'bmp'
   other options are remaining.

2. make shell script
   make batch scipt to run simpleISP with each sub-directory's param.json 


