import os
import argparse
def main(args):
   path=os.path.join(os.path.abspath(os.path.dirname(__file__)), f'info_{args.halide}.txt')
   os.remove(path)
if __name__=='__main__':

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
   
    parser.add_argument('-halide', type=str, default='F', help='Type of halide')
    args = parser.parse_args()
    main(args)
