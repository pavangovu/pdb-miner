from biopandas.pdb import PandasPdb
import re
import os
import argparse

def main(args):
    def filter(): 
            with open(f'info_{args.halide}.txt', 'a') as w:
                 w.write(f'{header}!{resolution}!{model_name}\n')

            with open(f'info_{args.halide}.txt', 'r') as f:
              global base_list
              base_list=[]
              text=f.readlines()
              for line in text:
                    base=list(line.split('!'))
                    if base[0] in base_list:
                               index=base_list.index(base[0])+1
                               if base[1]<base_list[index]:
                                   base_list[index-1]=base[0]
                                   base_list[index]=base[1]
                                   base_list[index+1]=base[2]
                    else:
                            base_list+=base
              print(base_list)

            with open(f'filtered_{args.halide}.txt', 'w') as f:
        
                                print(*base_list, file=f, sep="\n")

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
   
    try:
         header=re.search("COMPND\s+2\s+MOLECULE\:\s+(.+)\S+", struct.pdb_text).group(1)
    except:
         header=model_name
        
    permition=filter()

if __name__=='__main__':

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('-input', type=str,
#                         help='Path to input file.')
                        help='PDB id or PDB structure in .ent format.')
    parser.add_argument('-input_type', type=str, help='Pass your input type.')
    parser.add_argument('-halide', type=str, default='F', help='Type of halide')
    args = parser.parse_args()
    main(args)

