from biopandas.pdb import PandasPdb
import re
import os
import argparse
import sys


def filter(model_name,resolution,header,experiment): 

  with open(args.output_info, 'a') as w:
       w.write(f'{header}!{resolution}!{model_name}!{experiment}\n')

  with open(args.output_info, 'r') as f:
    global base_list
    base_list=[]
    text=f.readlines()
    number_inputs_before_filter=len(text)
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
    # print(base_list)

  with open (args.output_filtred, 'w') as f:

                      # print(*base_list, file=f, sep="\n")
                      line = '\n'.join(base_list)
                      f.write(line)
                      number_x_ray=base_list.count('X-RAY DIFFRACTION\n')
                      number_NMR=base_list.count('NMR\n')
                      number_inputs_after_filter=int(len(base_list)/4)
                      low_resolution=0
                      for i in range(1,len(base_list),4):
                         if float(base_list[i])>2:
                             low_resolution+=1
  return([number_inputs_before_filter,number_inputs_after_filter,number_x_ray,number_NMR,low_resolution])

def main(args):

   
    if args.input_type == 'pdb_id':

      try:
        struct = PandasPdb().fetch_pdb(args.input)
        model_name = args.input
      except:
        sys.exit()
        
    elif args.input_type == 'structure':

        struct = PandasPdb()
        struct = struct.read_pdb(args.input)
        model_name = re.search('[\d\w]+$', struct.header).group()
        
    
    try:
          resolution = float(re.search("REMARK\s+2\s+RESOLUTION\.\s+(\d+\.\d+)", struct.pdb_text).group(1))
    except:
          resolution = 100
    
    
    try:
         header=re.search("COMPND\s+2\s+MOLECULE\:\s+(.+)\S+", struct.pdb_text).group(1)
    except:
         header=model_name
    try:
       experiment =re.search("REMARK\s+\d+\s+EXPERIMENT TYPE\s+\:\s+(.+\w)\s+", struct.pdb_text).group(1)
    except:
       experiment='NMR'
    print(f'{args.input}: {experiment}')
    permition=filter(model_name,resolution,header,experiment)
    print(f'entries have been viewd: {permition[0]}, entries have been selected: {permition[1]}, entries with X-RAY: {permition[2]}, entries with NMR: {permition[3]}, entries with resolution > 2: {permition[4]}')
if __name__=='__main__':

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('-input', type=str,
#                         help='Path to input file.')
                        help='PDB id or PDB structure in .ent format.')
    parser.add_argument('-input_type', type=str, help='Pass your input type.')
    # parser.add_argument('-halide', type=str, default='F', help='Type of halide')
    parser.add_argument('-output_info', type=str)
    parser.add_argument('-output_filtred', type=str)

    args = parser.parse_args()
    main(args)

