##


def gen(wildcards):
	return [f'data/context/{i}_context' for i in ['BR', 'CL', 'F', 'I']]
	# return [f'data/context/{i}_context/{i}_{j}' for i in ['BR', 'CL', 'F', 'I'] for j in ['HIGH', 'MODERATE', 'LOW']]

rule all:
	input: gen
	output: touch('.status')

rule fetch_structures:
	input: "data/pdb_ID/pdb_entries_{halide}.txt"
	output: directory('data/structures/{halide}_struct')
	shell: 
		# 'mkdir data/strucures/{halide} |'
		'while read i; do python3.6 scripts/fetch_pdb.py -pdb_id $i -output {output} -format pdb; done < {input}'

rule filter:
	input: "data/structures/{halide}_struct"
	output: 
		info="data/info/info_{halide}.txt",
		filtred="data/filtered_pdb_ID/filtred_{halide}.txt"
	shell: "for i in $(ls {input}); do python scripts/filter.py -input {input}/$i -input_type structure -output_info {output.info} -output_filtred {output.filtred}; done"


rule get_context_data:
	input: 
		struct='data/structures/{halide}_struct/',
		filter='data/filtered_pdb_ID/filtred_{halide}.txt'
	output: directory('data/context/{halide}_context')
	shell: "for i in $(ls {input.struct}); do python3.6 scripts/get_context.py -input_struct $i -input_filter {input.filter} -input_type structure -angstrem_radius 5 -output {output} -C 2; done"

