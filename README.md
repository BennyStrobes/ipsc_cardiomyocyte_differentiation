# iPSC cardiomyocyte differentiation

High-resolution analysis of temporal dynamics of genetic effects on gene expression throughout cardiomyocyte differentiation.

Contact Ben Strober (bstrober3@gmail.com) with any questions.

## Running the code

The analysis/code is divided into three sections (each of which is a subdirectory):
1. 'preprocess_expression': Aligns reads and performs expression quantification.
2. 'standard_eqtl_calling': Runs the WASP Combined Haplotype Test (CHT) in each time step, independently.
3. 'dynamic_eqtl_calling': Identifies linear and non-linear dynamic eQTLs.

Each section depends on the previous, and therefore must be run sequentially, in series.


## Computer cluster

This pipeline was designed to run on midway2 rcc.

## Authors

* **Ben Strober** -- [BennyStrobes](https://github.com/BennyStrobes)
