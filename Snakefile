

rule all:
	input: 'plots/'
	output: touch('.status')


rule fetch_structures:
	input: "data/pdb_ID/pdb_entries_{halide}.txt"
	output: directory('data/structures/{halide}_struct')
	priority: 100
	shell: 
		'''
		while read i; do python scripts/fetch_structures.py -pdb_id $i -output {output} &>/dev/null; done < {input}
		i=$(echo {input})
		echo "=============="
		echo "$i: $(cat {input} | wc -l) structures were obtained."
		echo "=============="
		'''

rule filter:
	input: "data/structures/{halide}_struct"
	output: 
		info="data/info/info_{halide}.txt",
		filter="data/filtered_pdb_ID/filtered_{halide}.txt"
	priority: 10
	threads: 4
	shell: "for i in $(ls {input}); do python scripts/filter.py -input {input}/$i -input_type structure -output_info {output.info} -output_filtred {output.filter}; done"


rule get_context_data:
	input: 
		struct='data/structures/{halide}_struct',
		filter='data/filtered_pdb_ID/filtered_{halide}.txt'
	output: directory('data/context/{halide}_context')
	threads: 4
	shell: "python scripts/get_context.py -input_struct {input.struct} -input_filter {input.filter} -input_type structure -angstrem_radius 5 -output {output} -C 2 -input_ligands y -prot y"


rule combine_final_data:
	input:	expand("data/context/{halide}_context", halide=['BR', 'CL', 'F', 'I'])
	output: "full_data.tsv"
	shell: 'python scripts/parse_context.py -f {input} -input data/context/\*/\* -output {output}'


rule figures_plotting:
	input: "full_data.tsv"
	output: directory('plots/')
	shell: 'Rscript --vanilla scripts/figs.R {input} {output} &>/dev/null'


