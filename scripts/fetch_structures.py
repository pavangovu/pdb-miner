from Bio.PDB import *
import argparse
import os
import glob
import sys


def main(args):
    
    try:
        os.mkdir(args.output)
    except:
        pass
    
    try:
        pdbl = PDBList()
        pdbl.retrieve_pdb_file(args.pdb_id, pdir=args.output, file_format=args.format)
    except:
        print(f'Downloading error. {args.pdb_id} is not in the PDB.')
        pass

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('-pdb_id', type=str,
                        help='Fetch PDB structure file from PDB server, and store it locally.')
    parser.add_argument('-output', type=str,
                    help='Path to output directory.')
    parser.add_argument('-format', type=str, default='pdb', help='Output file format.')
    args = parser.parse_args()

    if f'pdb{args.pdb_id.lower()}' in glob.glob(args.output.lower().split('.')[0]):
        print(f'{args.pdb_id} is already in your local directory.')
        sys.exit()

    main(args)
