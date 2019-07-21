##


def gen(wildcards):
	return [f'data/strucures/{i}_struct/' for i in ['BR', 'CL', 'F', 'I']]

rule all:
	input: gen
	output: touch('.status')

rule filter:
	input: "data/pdb_ID/pdb_entries_{halide}.txt"
	output: "data/info_{halide}.txt"
	shell: "while read i; do python3.6 scripts/filter.py -input $i -input_type pdb_id -output_info {output}; done < {input}"

rule fetch_structures:
	input: "data/info_{halide}.txt"
	output: directory('data/strucures/{halide}_struct/')
	shell: 
		# 'mkdir data/strucures/{halide} |'
		'while read i; do python3.6 scripts/fetch_pdb.py -pdb_id $i -output {output} -format pdb; done < {input}'

# rule get_context_data:
# 	input: 'data/strucures/{halide}_struct/'
# 	output:
# 	shell:

