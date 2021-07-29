#!/usr/bin/python3

# USAGE: python nameOfThisScript.py
# This script takes some user inputs (i.e. blast results, Zotu/otu table)
# and thresholds defined by user (qCov, %identity, difference between %identity)
# Then it will link it to an NCBI taxonomy file and Zotu/Otu table file

# By Mahsa Mousavi-Derazmahalleh ; Python V3

import subprocess
import sys
# import supporting functions used in the code below
import working_function as wf

# define commands for downloading taxonomy database from NCBI and saving it in a folder with the current date
# if folder exists, it overwrites it
getLineage = '''
d="$(date +"%d-%m-%Y")" 
mkdir ${d}_taxdump
cd ${d}_taxdump 
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip 
unzip new_taxdump.zip 
mv rankedlineage.dmp ..  
cd .. 
tr -d "\t" < rankedlineage.dmp > rankedlineage_tabRemoved.dmp 
rm -f ${d}_taxdump/*
rm -f rankedlineage.dmp
'''

# create a child bash process
process = subprocess.Popen('/bin/bash', stdin=subprocess.PIPE,
                           stdout=subprocess.PIPE, universal_newlines=True)
# execute getLineage command and print output on the screen
out, err = process.communicate(getLineage)
out


# setting up input parameters
tableFile = str(sys.argv[1])
blastFile = str(sys.argv[2])
qCovLim = sys.argv[3]
pidLim = sys.argv[4]
pid_diffCut = sys.argv[5]

# define output file, open it and redirect standard output to this file
def_output = sys.stdout
output = str(sys.argv[6])
f = open(output, 'w')
sys.stdout = f

# define main function of the script


def main():

    # link filtered blast results with taxonomy information
    # results are saved in a temporary file called 'interMediate_res.tab'
    wf.link_TaxFilblast(blastFile, pid_diffCut, qCovLim, pidLim)

    # filter blast result file to get only unique hits based on taxonomy id,
    # and filter them based on query coverage, percentage identity and difference of them (Diff)
    n = wf.filterBlast(blastFile, pid_diffCut, qCovLim, pidLim)

    # link zotu table and store the values (table) and header (lable)
    table = wf.link_OTUtable(tableFile)['tableDict']
    lable = wf.link_OTUtable(tableFile)['lableDict']

    # print the header to the output file
    for lab in lable:
        print("domain\tphylum\tclass\torder\tfamily\tgenus\tspecies\tOTU\tnumberOfUnq_BlastHits" +
              '\t' + str('\t'.join(lable[lab])))

    # open filtered blast results temporary file
    with open("interMediate_res.tab", "r") as file:

        acc = [[]]  # make an empty list to store data
        i = 0       # define flag used to detect first line of the temporary file
        linenr = 0  # current line number
        loccnt = 6  # local column counter

        for line in file:
            # This will remove not helpfull assignments as set inside '' in below line,
            # but will still keep them in intermediate file for record
            if not any(value in line for value in ('uncultured bacterium', 'uncultured organism', 'uncultured prokaryote', 'uncultured eukaryote', 'uncultured marine eukaryote', 'uncultured marine picoeukaryote', 'uncultured Chlorophyta', 'marine metazoan environmental sample', 'zooplankton environmental sample', 'uncultured fungus', 'uncultured metazoan', 'eukaryote marine clone', 'eukaryotic picoplankton environmental sample', 'uncultured alveolate', 'uncultured ciliate', 'uncultured Banisveld eukaryote', 'uncultured chlorophyta', 'uncultured prasinophyte', 'uncultured Rhizaria', 'uncultured stramenopile', 'uncultured microeukaryote', 'uncultured marine dinoflagellate', 'uncultured prasinophyte', 'uncultured Cryptophyta', 'Stramenopiles sp. MAST', 'eukaryote clone OLI11008', 'marine Metazoa environmental sample', 'Dinophyta sp.', 'invertebrate environmental sample', 'uncultured eukaryote')):
                # read line to the 'll' list
                ll = line.strip().split('\t')
                # for the first line add value of the 'll' to the 'acc' list
                if i == 0:
                    acc[0] = ll
                    i = 1  # set flag to 1 (processed first line)

                # for the next lines
                else:
                    # if the column 'OTU' is equal to the 'OTU' of the last line in acc
                    # compare remaining columns (number 0 - 6)
                    if ll[7] == acc[linenr][7]:
                        loccnt = 6  # reset column counter to the value 6

                        # for each column
                        while loccnt > 0:
                            # if the values are not equal, or one of them is empty,
                            # replace the column value in the last line of acc with phase "dropped"
                            # and repeat replacement for all columns after investigated column
                            # (current column to column number 6)
                            if ll[loccnt] != acc[linenr][loccnt] and (ll[loccnt] != '' or acc[linenr][loccnt] != ''):
                                c = 6
                                while c >= loccnt:
                                    acc[linenr][c] = "dropped"
                                    c -= 1

                            # if values are equal and not empty
                            # stop comparison process for this column
                            elif ll[loccnt] == acc[linenr][loccnt] and ll[loccnt] != '' and acc[linenr][loccnt] != '':
                                break

                            # go to previous column
                            loccnt -= 1

                    # if 'OTU' is not equal, append it to the 'acc' list and move to next line
                    else:
                        acc.append(ll[:])
                        linenr += 1

    # linking of the Zotu with filtering and LCA taxonomy results
    for x in acc:
        # If the Zotu was found also in the dictionary of the table
        if x[7] in table:
            # print it to the results file
            print('\t'.join(x[:8]) + '\t' + str(n['vals']
                                                [x[7]]) + '\t' + '\t'.join(table[x[7]]))

# end of the main function of the script


# execute main function of the script
main()

# close output file
f.close()

# restore standard output to default
sys.stdout = def_output

# cleaning up temporary files
cmd = '''
d="$(date +"%d-%m-%Y")"
mv rankedlineage_tabRemoved.dmp ${d}_taxdump
'''
process = subprocess.Popen('/bin/bash', stdin=subprocess.PIPE,
                           stdout=subprocess.PIPE, universal_newlines=True)
out, err = process.communicate(cmd)
out
