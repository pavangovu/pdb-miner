# Halide sites  
## Dependencies
You should install it manually: 
* [Anaconda](https://www.digitalocean.com/community/tutorials/how-to-install-anaconda-on-ubuntu-18-04-quickstart)
* [Freesasa](https://freesasa.github.io/)  
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
1. **fetch_structures.py** - obtain structures from PDB;  
2. **filter.py** - applies filters to PDB structures;  
3. **get_context.py** - applies filters to halide binding sites;   
4. **parse_context.py** - aggregates data from all halides, then return 2 output files:
  * aggregated information about (model_name, distances, angels, fASA, atom_type, residue_type, etc.);
  * aggregated information about compositions of amino acids.  
5. **figs.R** - makes graphic report.  
