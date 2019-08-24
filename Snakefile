##


# def gen(wildcards):
# 	# return [f'data/context/{i}_context' for i in ['BR', 'CL', 'F', 'I']]
# 	return [f'data/context/{i}_context/{i}_{j}' for i in ['BR', 'CL', 'F', 'I'] for j in ['HIGH', 'MODERATE', 'LOW']]

import os


# try:
# 	[os.makedirs(f'data/structures/{i}_struct') for i in ['BR', 'CL', 'F', 'I']]
# except:
# 	pass

rule all:
	input: "full_data.tsv"
	output: touch('.status')


rule fetch_structures:
	input: "data/pdb_ID/pdb_entries_{halide}.txt"
	output: protected(directory('data/structures/{halide}_struct'))
	priority: 100
	shell: 
		'while read i; do python scripts/fetch_pdb.py -pdb_id $i -output {output} -format pdb &>/dev/null; done < {input} | '
		'echo "==============" |'
		'echo "Structures were obtained.'

rule filter:
	input: "data/structures/{halide}_struct"
	output: 
		info="data/info/info_{halide}.txt",
		filter="data/filtered_pdb_ID/filtered_{halide}.txt"
	priority: 10
	threads: 4

	shell: "for i in $(ls {input}); do python3.6 scripts/filter.py -input {input}/$i -input_type structure -output_info {output.info} -output_filtred {output.filter}; done"

rule get_context_data:
	input: 
		struct='data/structures/{halide}_struct',
		filter='data/filtered_pdb_ID/filtered_{halide}.txt'
	output: directory('data/context/{halide}_context')
	threads: 4
	shell: "python3.6 scripts/get_context.py -input_struct {input.struct} -input_filter {input.filter} -input_type structure -angstrem_radius 5 -output {output} -C 2"



rule combine_final_data:
	input:	expand("data/context/{halide}_context", halide=['BR', 'CL', 'F', 'I'])
	output: "full_data.tsv"
	shell: 'python scripts/parse_context.py -f {input} -input data/context/\*/\* -output {output}'

