from biopandas.pdb import PandasPdb
import re
import os
import argparse


def main(args):

    if args.input_type == 'pdb_id':

        struct = PandasPdb().fetch_pdb(args.input)
        model_name = args.input

    elif args.input_type == 'structure':

        struct = PandasPdb()
        struct = struct.read_pdb(args.input)
        model_name = re.search('[\d\w]+$', struct.header).group()

    try:
        resolution = float(re.search("REMARK\s+2\s+RESOLUTION\.\s+([\d\.]+)\s\w+", struct.pdb_text).group(1))
    except:
        resolution = 100

    halide_atoms = struct.df['HETATM'][struct.df['HETATM']['atom_name'] == args.halide]

    dict_of_subsets = {}
    for i in halide_atoms.values:
        dist = struct.distance(xyz=tuple(i[11:14]), records=('ATOM'))
        atom_subset = struct.df['ATOM'][dist < args.angstrem_radius]
        dict_of_subsets[f'{model_name}:{resolution}:{i[3]}:{i[11]}:{i[12]}:{i[13]}'] =\
                        [(f'{j[3]}:{j[5]}:{j[11]}:{j[12]}:{j[13]}') for j in atom_subset.values]

    try:
        os.makedirs(args.output_dir)
    except:
        pass

    def write_output(sfx):

        with open(f'{args.output_dir}/{args.output_file_name}_{sfx}.tsv', 'a') as w:
            for k,v in dict_of_subsets.items():
                w.write(f'{k}\t')
                for i in range(len(v)):
                    if i == len(v)-1:
                        w.write(f'{v[i]}')
                    else:
                        w.write(f'{v[i]},')
                w.write('\n')

    if resolution <= 1.5:
        write_output('HIGH')
    elif resolution > 1.5 and resolution < 2.5:
        write_output('MODERATE')
    else:
        write_output('LOW')

if __name__=='__main__':

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('-input', type=str,
#                         help='Path to input file.')
                        help='PDB id or PDB structure in .ent format.')
    parser.add_argument('-input_type', type=str, help='Pass your input type.')
    parser.add_argument('-halide', type=str, default='CL', help='Type of halide')
    parser.add_argument('-angstrem_radius', type=int, default=5, help='Threshold radius in Å.')
    parser.add_argument('-output_file_name', type=str, 
                        help='Name of output file (root; suffixes will be put themselves).')
    parser.add_argument('-output_dir', type=str, 
                    help='Name of output dir.')
    args = parser.parse_args()

    main(args)
