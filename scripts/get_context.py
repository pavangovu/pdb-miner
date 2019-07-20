from biopandas.pdb import PandasPdb
import re
import os
import argparse
import numpy as np
import pandas as pd
import copy
import itertools
import freesasa


def ang(Coordinate): # get angles
             
                                Angles=[]
             
                                for  k in range(2,len(Coordinate)): 

                                     a=np.array(Coordinate[0])
                                     b=np.array(Coordinate[1])
                                     c=np.array(Coordinate[k])
                                     ab=b-a
                                     ac=c-a
                                     cosine_angle = np.dot(ab,ac)/(np.linalg.norm(ab)*np.linalg.norm(ac))
                                     cosine_angle=round(cosine_angle,15) # remove the error
                                     angle=np.arccos(cosine_angle) 
                                     if  np.degrees(angle) !=0:
                                         Angles+=[np.degrees(angle)] # fill the list of angles
                                     else:
                                         Angles+=[np.nan]
            
                                return(Angles)




def main(args):

      try:
        os.makedirs(args.output_dir)
      except:
        pass

      with open(f'../data/filtered_pdb_ID/filtered_{args.halide}.txt', 'r') as f:
          text=f.readlines()
          base_list=[]
          base_list+=text
          base_list=[line.rstrip() for line in base_list]
          base_list=[i for i in base_list if i !='']
          print (base_list)
          for i in range(2,len(base_list),3): 
              
              if args.input_type == 'pdb_id':

                   struct = PandasPdb().fetch_pdb(f'{base_list[i]}')
                   model_name = args.input
                   structure= freesasa.Structure(f'{base_list[i]}')

              elif args.input_type == 'structure':
                   struct = PandasPdb()
                   struct = struct.read_pdb(f'{args.input}/pdb{base_list[i].lower()}.ent')
                   print(f'{args.input}/pdb{base_list[i].lower()}.ent')
                   model_name = re.search('[\d\w]+$', struct.header).group()
                   structure= freesasa.Structure(f'{args.input}/pdb{base_list[i].lower()}.ent')
              try:
                        resolution = float(re.search("REMARK\s+2\s+RESOLUTION\.\s+([\d\(.)]+)\s\w+", struct.pdb_text).group(1))
              except:
                        resolution = 100
              print (resolution)
      
              result= freesasa.calc(structure)
              classArea =freesasa.classifyResults(result,structure)
              print(result.totalArea())

              halide_atoms = struct.df['HETATM'][struct.df['HETATM']['atom_name'] == args.halide]
              modern_df=struct.df['ATOM'] # make the subset 
              dict_of_subsets = {}
              for i in halide_atoms.values:
                          global Coordinate
                          Coordinate=[] # list of coordinates М[0]= halide coordinates М[1] the nearest atom's coordinates
                          dist = struct.distance(xyz=tuple(i[11:14]), records=('ATOM'))
                          modern_df['dist']=dist # add distanse to subset
                          

                          if args.C == 1:

                                 modern_subset =modern_df[modern_df.dist < args.angstrem_radius] # halide neighbors
                      
                          
                          elif args.C == 2:

                                 modern_df1=modern_df[modern_df.dist < args.angstrem_radius]
                                 modern_subset=modern_df1.loc[~modern_df1['element_symbol'].isin(['C','H'])]

        
                          elif args.C == 3:
                                 
                                 modern_df1 =modern_df[modern_df.dist < args.angstrem_radius]
                                 modern_subset=modern_df1.loc[~modern_df1['element_symbol'].isin(['C'])]

                          elif args.C == 4:
                                 modern_df1 =modern_df[modern_df.dist < args.angstrem_radius]
                                 atoms=list(map("".join,itertools.permutations('BDEGHZ',1)))
                                 atoms1=list(map("".join,itertools.permutations('BDEGHZ123',2)))
                                 atoms2=atoms+atoms1
                                 atom2=['C'+ atom for atom in atoms2] +['O','C']
                                 modern_subset=modern_df1.loc[~modern_df1['atom_name'].isin(atom2)]

                                
                          
                          xyz=i[11:14] # halide coordinate
                          Coordinate+=[xyz] # add coordinates
                          try:
                                    nearest=modern_subset.loc[modern_subset['dist']==min(modern_subset['dist'])].values[0][[3,5,11,12,13,20,21]] #define the nearest atom
            
                                    Coordinate+=[nearest[2:5]] # add the nearest atom's coordinates
                          except:
                               continue
                          for n in modern_subset.values:
                                 xyz2= n[11:14]
                                 Coordinate+=[xyz2] # add coordinates
        
                          modern_subset1=modern_subset.copy(deep=True)
                          modern_subset1['angles']=ang(Coordinate) # add angles to subset
                          #modern_subset1=modern_subset1.loc[modern_subset1['angles'] != 0] # delete rows  with angles=0
                          #{nearest[0]}:{nearest[1]}:{"%.3f"% nearest[6]}:{np.nan}
                          dict_of_subsets[f'{model_name}:{"%.3f"% result.totalArea()}:{resolution}:{i[3]}:{nearest[5]}'] =\
                          [(f'{j[3]}:{j[5]}:{"%.3f"% j[21]}:{"%.3f"% j[22]}') for j in modern_subset1.values]
 

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
      path=os.path.join(os.path.abspath(os.path.dirname(__file__)), 'pdb_one_halide.txt')
      os.remove(path)

if __name__=='__main__':

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('-input', type=str,
#                         help='Path to input file.')
                        help='PDB id or PDB structure in .ent format.')
    parser.add_argument('-input_type', type=str, help='Pass your input type.')
    parser.add_argument('-halide', type=str, default='F', help='Type of halide')
    parser.add_argument('-angstrem_radius', type=int, default=5, help='Threshold radius in Å.')
    parser.add_argument('-output_file_name', type=str, 
                        help='Name of output file (root; suffixes will be put themselves).')
    parser.add_argument('-output_dir', type=str, 
                    help='Name of output dir.')
    parser.add_argument('-C', type=int, help='1-all atoms; 2-no C,H atoms; 3 - no C; 4 - no C=0 no C except CA)
    args = parser.parse_args()

    main(args)
