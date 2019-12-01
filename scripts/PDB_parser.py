import requests
import argparse

def main(args):
    symbols ={'F':'fluoride','CL':'chloride','BR':'bromide','I': 'iodide' }
    url = 'http://www.rcsb.org/pdb/rest/search'
    query_text = f'<?xml version="1.0" encoding="UTF-8"?>\
    <orgPdbQuery>\
    <version>B0907</version>\
    <queryType>org.pdb.query.simple.ChemCompNameQuery</queryType>\
    <description>Experimental Method Search: Experimental Method=SOLID-STATE NMR</description>\
    <comparator>Equals</comparator>\
    <name>{symbols[args.ian]}</name>\
    <polymericType>Any</polymericType>\
    </orgPdbQuery>'

    print("Querying RCSB PDB REST API...")

    header = {'Content-Type': 'application/x-www-form-urlencoded'}
  
    response = requests.post(url, data=query_text, headers=header)
    if response.status_code == 200:
        print(f'Found {int(len(response.text)/5)} PDB entries matching query.')
        with open(f'{args.output}', 'w') as w:
  
             w.write(response.text)
    else:
        print("Failed to retrieve results")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()    
    parser.add_argument('-ian', type=str, help = 'Enter list of atom names, such as fluoride...')
    parser.add_argument('-output', type=str, help = 'Output files')
    args = parser.parse_args()

    main(args)
