referencing 'ctest/parameters.txt'
referencing 'Test/task1.md'

### Making parameters 
1. making parameters.txt files for each subfolder of Test
2. change fields according to rules like below
BayerImg field: selecting one file path out of ctest/bayerimg
width and height: extract information from BayerImg file name
bpp: 16
bayer_img_csi2_style: false

bitstream_filename:  each folder path + its folder name + '.bin'
codec_report_filename:  each folder path + its folder name + '.rpt'
dec_cfa_filename = each folder path + 'bayer_dec'+each folder name + '.raw'

dwt_lv, ngrp, sgpcd, nlt, sel_ycc, Qp, Gb[0] ~ Gb[4]: extract values from folder names.
others: copy parameters from 'ctest/parametes.txt'.

### Making shell script
1. make a shell script named as run_batch_codec.sh
2. It executes 'ctest/build/CfaCodec' with each parameters in each subfolder.
 
