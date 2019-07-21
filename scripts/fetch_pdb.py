from Bio.PDB import *
import argparse
import os

def main(args):
    
    try:
        os.mkdir(args.output)
    except:
        pass
    
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(args.pdb_id, pdir=args.output, file_format=args.format)

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('-pdb_id', type=str,
                        help='Fetch PDB structure file from PDB server, and store it locally.')
    parser.add_argument('-output', type=str,
                    help='Path to output directory.')
    parser.add_argument('-format', type=str, default='pdb', help='Output file format.')
    args = parser.parse_args()
    main(args)
