# ipsc preprocess pipeline

This pipeline processes/quantifies the fastq files to organized allelic count matrices & quantile normalized expression matrices. 

## Running the code

Preprocessing can be run within 'preprocess_driver.key.sh'. 
'preprocess_driver.key.sh' splits pipeline into 6 steps that can be run in series (all the user has to do is remove the if false; then code arround the desired part and submit the current job). 
Comments within 'preprocess_driver_key'.sh' explain what each section does. But briefly:


## Computer cluster

This pipeline was written to run on midway rcc

## Authors

* **Ben Strober** -- [BennyStrobes](https://github.com/BennyStrobes)
