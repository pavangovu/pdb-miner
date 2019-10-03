# Halide sites  
## Dependencies
You should install it manually: 
* [Anaconda](https://www.digitalocean.com/community/tutorials/how-to-install-anaconda-on-ubuntu-18-04-quickstart)
* [FreeSASA](https://freesasa.github.io/)  
## Preparation  
Create your conda environment from our smk.yml  
```
conda env create -f smk.yml --name smk
```
Activate your environment (some of this should work)  
```
conda activate smk
```
or
```
source activate smk
```

Make sure everything works on `--dryrun`.
```
snakemake all -n
```
## How to use
```
snakemake all
```
## DAG of jobs  
![alt text](dag.svg)  
## Scripts description  
1. rule fetch_structures (**fetch_structures.py**) - obtain structures from PDB;  
2. rule filter (**filter.py**) - applies filters to PDB structures;  
3. rule get_context_data (**get_context.py**) - applies filters to halide binding sites;   
4. rule combine_final_data (**parse_context.py**) - aggregates data from all halides, then return 2 output files:
  * aggregated information about (model_name, distances, angels, fASA, atom_type, residue_type, etc.);
  * aggregated information about compositions of amino acids.  
5. rule figures_plotting (**figs.R**) - makes graphic report.  
## References  
1. Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web. <https://anaconda.com>.  
2. 1. Mitternacht, S. FreeSASA: An open source C library for solvent accessible surface area calculations. F1000Research 5, 189 (2016).  
