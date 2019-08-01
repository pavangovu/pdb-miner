from biopandas.pdb import PandasPdb
import re
import os
import argparse
import numpy as np
import pandas as pd
import copy
import itertools
# import freesasa
import subprocess


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
 
def site_filter_ang(N,modern_subset1,sites_ang_sort):
  for k in range(0,N):
           same_degree=0
           for g in range(len(modern_subset1['angles'])):
                  if abs(sites_ang_sort.iloc[g][N]-sites_ang_sort.iloc[g][k])<4:
                          same_degree+=1
                  if (same_degree==len(modern_subset1['angles']) and (same_degree!=1)):
                          return True

def site_filter_dist(N,modern_subset1,sites_dist_sort):
   for k in range(0,N):
    same_degree=0
    for g in range(len(modern_subset1['dist'])):
      if abs(sites_dist_sort.iloc[g][N]-sites_dist_sort.iloc[g][k])<0.5:
        same_degree+=1
      if same_degree==len(modern_subset1['dist']):
        return True



def main(args):

  try:
    os.makedirs(args.output_dir)
  except:
    pass

  with open(args.input_filter, 'r') as f:
    text=f.readlines()
    base_list=[]
    base_list+=text
    base_list=[line.rstrip() for line in base_list]
    base_list=[i for i in base_list if i !='']
    # print (base_list)
    for i in range(2,len(base_list),3): 
        
      if args.input_type == 'pdb_id':

       struct = PandasPdb().fetch_pdb(f'{base_list[i]}')
       model_name = args.input
       with open('current_pdb.txt', 'w') as w:
                       w.write(f'{struct.pdb_text}')        
       model_name = base_list[i]
       # print(model_name)
       with open('current_pdb.txt', 'r') as w:
                       f=w.readlines()
       for k in range(len(f)):
             if (f[k][0:6]=='HETATM') and (f[k][16:20]==' HOH' or f[k][16:20]=='AHOH' or f[k][16:20]=='BHOH'):
                       f[k]=''
           

      elif args.input_type == 'structure':
        struct = PandasPdb()

        struct = struct.read_pdb(f'{args.input_struct}/pdb{base_list[i].lower()}.ent')
        # print(f'{args.input}/pdb{base_list[i].lower()}.ent')
        model_name = re.search('[\d\w]+$', struct.header).group()
        with open (f'{args.input_struct}/pdb{base_list[i].lower()}.ent', 'r') as pdb1:

          struct = struct.read_pdb(f'{args.input_struct}/pdb{base_list[i].lower()}.ent')
          # print(f'{args.input}/pdb{base_list[i].lower()}.ent')
          model_name = re.search('[\d\w]+$', struct.header).group()
        with open (f'{args.input_struct}/pdb{base_list[i].lower()}.ent', 'r') as pdb1:

                      f=pdb1.readlines()
        for k in range(len(f)):
              if (f[k][0:6]=='HETATM') and (f[k][16:20]==' HOH' or f[k][16:20]=='AHOH' or f[k][16:20]=='BHOH'):
                      f[k]=''
        try:
                  resolution = float(re.search("REMARK\s+2\s+RESOLUTION\.\s+(\d+\.\d+)", struct.pdb_text).group(1))
        except:
                  resolution = 100
        # print(resolution)

       

        halide_type = args.input_filter.split('/')[-1].split('_')[-1].split('.')[0]
        # print(halide_type)

        halide_atoms = struct.df['HETATM'][struct.df['HETATM']['atom_name'] == halide_type]
       
        modern_df=struct.df['ATOM'] # make the subset 
        dict_of_subsets = {}
        S=0 # halide counter for using freesasa
        sites_ang=pd.DataFrame({'number': [np.nan for k in range(0,100)]}).dropna(axis=1, how='all')
        sites_dist=pd.DataFrame({'number': [np.nan for k in range(0,100)]}).dropna(axis=1, how='all')
        N=-1 #sites counter for sites filter
        for i in halide_atoms.values:
                    
                    halide_atoms.index=np.arange(len(halide_atoms))
                    Halide_humber= halide_atoms[halide_atoms.index==S].values[0][1]
                    N+=1
                    S+=1
                    f1=copy.deepcopy(f)
                    for k in range(len(f1)):
                                if ((f1[k][0:6]=='HETATM') and (int(f1[k][6:11])!=int(Halide_humber))):
                                  if (f1[k][77:78]==halide_type) or (f1[k][76:78]==halide_type):
                                    f1[k]=''
                    f1=[x for x in f1 if x]
                    

                    # with open ('pdb_one_halide.txt', 'w') as pdb:

                    with open ('../pdb_one_halide.txt', 'w') as pdb:

                              #print(*f1,file=pdb,sep='\n')
                              line = '\n'.join(f1)
                              pdb.write(line)

                    out=subprocess.Popen(["freesasa", "../pdb_one_halide.txt", "-H", "--select", f'asa, symbol {halide_type}'],
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT,
                                       encoding='utf-8')
                    res=out.communicate()
                    if res[1]==None:
                         asa=re.search("asa\s+\:\s+(\d+\.\d+)", res[0]).group(1)
                         # print(asa)
                    else:
                         asa=np.nan
                         # print(asa)
                   
                   


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

                    modern_subset2=modern_subset.copy(deep=True)
                    modern_subset2.loc[:,['x_coord','y_coord','z_coord']]=modern_subset[['x_coord','y_coord','z_coord']]-[i[11],i[12],i[13]]      
                    
                    xyz=i[11:14] # halide coordinate
                    for k in range(len(xyz)):
                               xyz[k]=0
                 
                    Coordinate+=[xyz] # add coordinates
                    try:
                              nearest=modern_subset.loc[modern_subset['dist']==min(modern_subset['dist'])].values[0][[3,5,11,12,13,20,21]] #define the nearest atom
      
                              Coordinate+=[nearest[2:5]] # add the nearest atom's coordinates
                    except:
                         continue
                    for n in modern_subset.values:
                           xyz2= n[11:14]
                           Coordinate+=[xyz2] # add coordinates
  
<<<<<<< HEAD

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

                          modern_subset2=modern_subset.copy(deep=True)
                          modern_subset2.loc[:,['x_coord','y_coord','z_coord']]=modern_subset[['x_coord','y_coord','z_coord']]-[i[11],i[12],i[13]]      
                          
                          xyz=i[11:14] # halide coordinate
                          for k in range(len(xyz)):
                                     xyz[k]=0
                       
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

                          sites_ang[f'angles_{N}']=modern_subset1['angles']
                          sites_dist[f'dist_{N}']=modern_subset1['dist']
                         
                          sites_ang_sort=pd.DataFrame(np.sort(sites_ang, axis=0),columns=sites_ang.columns)
                          sites_dist_sort=pd.DataFrame(np.sort(sites_dist, axis=0),columns=sites_dist.columns)

                          if (site_filter_ang(N,modern_subset1['angles'],sites_ang_sort)==1) and (site_filter_dist(N,modern_subset1['dist'],sites_dist_sort)==1):

                               print(f'{halide_type} atom is skipped, similar haligen site have already been got')
                               continue



                          #modern_subset1=modern_subset1.loc[modern_subset1['angles'] != 0] # delete rows  with angles=0
                          #{nearest[0]}:{nearest[1]}:{"%.3f"% nearest[6]}:{np.nan}
                          dict_of_subsets[f'{model_name}:{asa}:{resolution}:{i[3]}:{nearest[5]}'] =\
                          [(f'{j[3]}:{j[5]}:{"%.3f"% j[21]}:{"%.3f"% j[22]}') for j in modern_subset1.values]
 

              def write_output(sfx):
                try:
                  os.makedirs(f'data/context/{halide_type}_context')
                except:
                  pass

                with open(f'{args.output}/{halide_type}_context_{sfx}.tsv', 'a') as w:
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
      # path=os.path.join(os.path.abspath(os.path.dirname(__file__)), '../pdb_one_halide.txt')
      # os.remove(path)

                    modern_subset1=modern_subset.copy(deep=True)
                    modern_subset1['angles']=ang(Coordinate) # add angles to subset

                    sites_ang[f'angles_{N}']=modern_subset1['angles']
                    sites_dist[f'dist_{N}']=modern_subset1['dist']
                   
                    sites_ang_sort=pd.DataFrame(np.sort(sites_ang, axis=0),columns=sites_ang.columns)
                    sites_dist_sort=pd.DataFrame(np.sort(sites_dist, axis=0),columns=sites_dist.columns)
                    if site_filter_ang(N,modern_subset1['angles'],sites_ang_sort)=1 and site_filter_dist(N,modern_subset1['dist'],sites_dist_sort)=1:

                         print(f'{halide_type} atom is skipped, similar haligen site have already been got')
                         continue



                    #modern_subset1=modern_subset1.loc[modern_subset1['angles'] != 0] # delete rows  with angles=0
                    #{nearest[0]}:{nearest[1]}:{"%.3f"% nearest[6]}:{np.nan}
                    dict_of_subsets[f'{model_name}:{asa}:{resolution}:{i[3]}:{nearest[5]}'] =\
                    [(f'{j[3]}:{j[5]}:{"%.3f"% j[21]}:{"%.3f"% j[22]}') for j in modern_subset1.values]


        def write_output(sfx):
          try:
            os.makedirs(f'data/context/{halide_type}_context')
          except:
            pass

          with open(f'{args.output}/{halide_type}_context_{sfx}.tsv', 'a') as w:
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
  # path=os.path.join(os.path.abspath(os.path.dirname(__file__)), '../pdb_one_halide.txt')
  # os.remove(path)


if __name__=='__main__':

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('-input_filter', type=str)
    parser.add_argument('-input_struct', type=str,
#                         help='Path to input file.')
                        help='PDB id or PDB structure in .ent format.')
    parser.add_argument('-input_type', type=str, help='Pass your input type.')
    # parser.add_argument('-halide', type=str, default='F', help='Type of halide')
    parser.add_argument('-angstrem_radius', type=int, default=5, help='Threshold radius in Å.')
    parser.add_argument('-output', type=str, 
                        help='Name of output file (root; suffixes will be put themselves).')
    # parser.add_argument('-output_dir', type=str, 
    #                 help='Name of output dir.')
    parser.add_argument('-C', type=int, help='1-all atoms; 2-no C,H atoms; 3 - no C; 4 - no C=0 no C except CA')
    args = parser.parse_args()

    main(args)
