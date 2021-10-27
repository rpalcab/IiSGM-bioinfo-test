import pandas as pd
import io
import os
import re

def read_vcf(path):
    ''' ( file_path ) --> Dataframe
    Introduce the path of the vcf file.'''

    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

def vcf_to_table(file_path):
    ''' ( file_path ) --> Dataframe
    Introduce the path of the vcf file.'''

    with open(file_path, 'r') as vcf_file:

        lines = []

        for line in vcf_file:
            if line.startswith('#CHROM'):
                lines.append(line)
            elif not line.startswith('#') and line.split('\t')[6] == 'PASS':
                lines.append(line)
                
        table = pd.read_csv(
                            io.StringIO(''.join(lines)),
                            dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                                'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT': str},
                            sep='\t'
                        ).rename(columns={'#CHROM': 'CHROM'})
        return table

def absence_presence_matrix(files_path):
    ''' ( directory_path ) --> Dataframe
    Introduce the path of the directory with the vcf files.'''

    # Create a empty dict
    dict_POS_REF_ALT = {}

    # Iterate through a given directories 
    for folderName, subfolders, filenames in os.walk(files_path):
        for filename in filenames:
            
            # Create the path of the vcf file and use the vcf_to_table function
            path = os.path.join(folderName, filename)
            table = vcf_to_table(path)
            
            # Create a list with the names of our samples, we will used it later tofill with ones or leave it empty
            list_samples = [filename.split('.')[0] for filename in filenames]

            # Create a dict with the name of all samples with value 0, to change 0 by 1 if that sample have that SNP
            list_samples_dict = dict(zip(list_samples, [0] * len(list_samples)))

            for row in range(len(table)):

                # We create the name of the SNP by its position, its reference base and its alternative base
                # We use the SNP name as the columns of the presence/absence matrix
                name_row =  str(table.POS[row]) + '_' + str(table.REF[row]) + '_' + str(table.ALT[row])

                # We create an exception that only occurs if the SNP have more than one alternative base
                if len(table.ALT[row].split(',')) > 1:
                    
                    # We iterate by the genotypes in GT (1/2), and if have an alternative base use it to create the name of the SNP
                    for gt in re.findall(r"[\w']+", table.iloc[:,-1][row].split(':')[0]):
                        if int(gt) > 0:
                            name_row = str(table.POS[row]) + '_' + str(table.REF[row]) + '_' + str(table.ALT[row].split(',')[int(gt) - 1])

                # If the name of the SNP does not exist in our dict create it and fill it with the dict with all the samples
                if name_row not in dict_POS_REF_ALT:
                    dict_POS_REF_ALT[name_row] = list_samples_dict

                # If a sample have this SNP change 0 by 1
                dict_POS_REF_ALT[name_row][filename.split('.')[0]] = 1

    # Turn our dict in to a table
    return pd.DataFrame(dict_POS_REF_ALT)

def pairwise_dist_table(files_path):
    ''' ( directory_path ) --> Dataframe
    Introduce the path of the directory with the vcf files.'''

    # Create a empty dict
    dict_POS_REF_ALT = {}
    dict_pairwise_dist = {}

    # Iterate through a given directories 
    for folderName, subfolders, filenames in os.walk(files_path):
        for filename in filenames:

            # Create the path of the vcf file and use the vcf_to_table function
            path = os.path.join(folderName, filename)
            table = vcf_to_table(path)

            dict_POS_REF_ALT[filename.split('.')[0]] = []

            for row in range(len(table)):
                
                # We create the name of the SNP by its position, its reference base and its alternative base
                name_row =  str(table.POS[row]) + '_' + str(table.REF[row]) + '_' + str(table.ALT[row])

                # We create an exception that only occurs if the SNP have more than one alternative base
                if len(table.ALT[row].split(',')) > 1:
                    
                    # We iterate by the genotypes in GT (1/2), and if have an alternative base use it to create the name of the SNP
                    for gt in re.findall(r"[\w']+", table.iloc[:,-1][row].split(':')[0]):
                        if int(gt) > 0:
                            name_row = str(table.POS[row]) + '_' + str(table.REF[row]) + '_' + str(table.ALT[row].split(',')[int(gt) - 1]) 

                dict_POS_REF_ALT[filename.split('.')[0]].append(name_row)

    for sample1 in dict_POS_REF_ALT:
        dict_pairwise_dist[sample1] = {}
        for sample2 in dict_POS_REF_ALT:
            dict_pairwise_dist[sample1][sample2] = abs(len(dict_POS_REF_ALT[sample1]) - len(dict_POS_REF_ALT[sample2]))

    return pd.DataFrame(dict_pairwise_dist)

def pairwise_dist_table_INDELS(files_path):
    ''' ( directory_path ) --> Dataframe
    Introduce the path of the directory with the vcf files.'''

    # Create a empty dict
    dict_POS_REF_ALT = {}
    dict_pairwise_dist = {}

    # Iterate through a given directories 
    for folderName, subfolders, filenames in os.walk(files_path):
        for filename in filenames:

            # Create the path of the vcf file and use the vcf_to_table function
            path = os.path.join(folderName, filename)
            table = vcf_to_table(path)

            dict_POS_REF_ALT[filename.split('.')[0]] = []

            for row in range(len(table)):

                # We create the name of the SNP by its position, its reference base and its alternative base
                name_row =  str(table.POS[row]) + '_' + str(table.REF[row]) + '_' + str(table.ALT[row])

                # We create an exception that only occurs if the SNP have more than one alternative base
                if len(table.ALT[row].split(',')) > 1:
                    
                    # We iterate by the genotypes in GT (1/2), and if have an alternative base use it to create the name of the SNP
                    for gt in re.findall(r"[\w']+", table.iloc[:,-1][row].split(':')[0]):
                        if int(gt) > 0:
                            name_row = str(table.POS[row]) + '_' + str(table.REF[row]) + '_' + str(table.ALT[row].split(',')[int(gt) - 1]) 

                # We create an exception for INDELS, when the alternative base be more than 1 of length create so many entries in the SNP dict as the length of the INDEL
                if len(name_row.split('_')[2]) > 1:
                    for number in range(len(str(table.ALT[row].split(',')[int(gt) - 1]))):
                        dict_POS_REF_ALT[filename.split('.')[0]].append(name_row + str(number))

                else:
                    dict_POS_REF_ALT[filename.split('.')[0]].append(name_row)

    # We compare a sample against another and count the result, that result will be the pairwise distance between the samples
    for sample1 in dict_POS_REF_ALT:
        dict_pairwise_dist[sample1] = {}
        for sample2 in dict_POS_REF_ALT:
            dict_pairwise_dist[sample1][sample2] = abs(len(dict_POS_REF_ALT[sample1]) - len(dict_POS_REF_ALT[sample2]))

    return pd.DataFrame(dict_pairwise_dist)