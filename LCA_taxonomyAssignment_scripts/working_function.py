#!/usr/bin/python3
# These are all functions necessary to run the LCA script
# By Mahsa Mousavi-Derazmahalleh ; Python V3

# This function is to split the lines in blast result file
def splitFile(ll):
    otuid = ll[0]       # otu id
    staxids = ll[2]     # taxonomy id
    pident = ll[6]      # percentage identity
    evalue = ll[18]     # evalue
    qcov = ll[20]       # query coverage
    sseqid = ll[1]      # subject sequence id
    commonames = ll[4]  # common name
    length = ll[7]      # length of alignment
    qlen = ll[8]        # query length
    slen = ll[9]        # sunject length
    return {'otuid': ll[0], 'staxids': staxids, 'pident': pident, 'evalue': evalue, 'qcov': qcov, 'sseqid': sseqid, 'commonames': commonames, 'length': length, 'qlen': qlen, 'slen': slen}


""" This is to filter blast result file to get only unique hits based on taxonomy id (so if there are few hits to the same taxonomy id, only the 1st one is kept) and then filter them more """


def filterBlast(filename, diff_lim, qCovThre, pidThre):
    taxidDict = {}       # a dictionary of otu as a key and taxonomy id as the value
    taxId_seen = set()   # keep taxonomy id which already seen
    unq = []             # for keeping the unique hits i.e. unique based on taxonomy id
    filunq = []          # an empty list to store result that pass initial qcov and percentage identity thresholds set by user
    vals = {}            # create a dictionary for otu>>>number of unique hits

    with open(filename, "r") as file:
        for line in file:
            bl = line.strip().split('\t')
            # store the result of the splitFile function in the value
            n = splitFile(bl)
            # To have multiple value for the same key
            taxidDict.setdefault(n['otuid'], []).append(n['staxids'])
            # make a combination for otu/taxid to later check their repetition
            notSeen = n['otuid'], n['staxids']

            # if the combination above was not seen before, append it to unq, and add that combination to taxId_seen set
            if notSeen not in taxId_seen:
                unq.append(bl[:])
                taxId_seen.add(notSeen)

    # making a dictionary to hold the number of unique hits per otu, where key is otu and value is number of unique hit for that
    for otid, taid in list(taxidDict.items()):
        vals[otid] = len(set(taid))

    # if qcov >= user specified threshold for qCovThre and percentage identity >= user specified threshold for pidThre append to filunq
    for i in unq:
        if (round(float(i[20]), 2) >= round(float(qCovThre), 2)) and (round(float(i[6]), 2) >= round(float(pidThre), 2)):
            filunq.append(i[:])

    c = 0  # current line in unq
    d = 1  # flag if unq was changed

    while d > 0:
        c = 0
        d = 0
        # This is to compare the first line with the next, and so on
        for x, y in zip(filunq, filunq[1:]):
            otuId_x = x[0]
            otuId_y = y[0]
            pIdent_x = float(x[6])
            pIdent_y = float(y[6])
            qCov_x = float(x[20])
            qCov_y = float(y[20])

            # if otu ids in the two lines are the same
            if otuId_x == otuId_y:
                # if query coverage in line1 and line2 were equal too
                if qCov_x == qCov_y:
                    # if absolute value for the difference between %identity of line1 vs line2 is > than the threshold set by user
                    if "{0:.3f}".format(abs(pIdent_x - pIdent_y)) > "{0:.3f}".format(float(diff_lim)):
                        # % identity of line one is bigger than line two
                        if pIdent_x > pIdent_y:
                            # remove the second line from unq (only keep the best between the two)
                            del filunq[c+1]
                            d = d + 1
                            break
                        else:
                            # remove the first line from unq (only keep the best between the two)
                            del filunq[c]
                            d = d + 1
                            break

                elif qCov_x > qCov_y:
                    del filunq[c+1]
                    d = d + 1
                    break

                elif qCov_x < qCov_y:
                    del filunq[c]
                    d = d + 1
                    break
            c = c + 1

    return {'unq': filunq, 'vals': vals}


# This function will make a dictionary for taxonomy information of ncbi
def taxonomy_dictionary(ncbi_taxonomy):

    taxonomy_Dict = {}

    # This will parse the ncbi taxonomy file to get the column of interest
    with open(ncbi_taxonomy, "r") as file:
        for line in file:
            l3 = line.strip().split('|')
            tax_id = l3[0]
            domain = l3[9]
            phylum = l3[7]
            classif = l3[6]
            order = l3[5]
            family = l3[4]
            genus = l3[3]
            species = l3[1]

            if "tax_id" not in line:
                taxonomy_Dict[tax_id] = [domain, phylum,
                                         classif, order, family, genus, species]

    return taxonomy_Dict


# This function will link filtered blast results with taxonomy information
# it takes four arguments including the filtering thresholds set by users, and will write the result to an intermediate file
def link_TaxFilblast(filename, diff_lim, qCovThre, pidThre):
    f = open('interMediate_res.tab', 'w')

    # holds the value from filterBlast function. which is a nested list
    filBlast = filterBlast(filename, diff_lim, qCovThre, pidThre)['unq']
    # a dictionary of taxonomyid and the value is all information available for it
    taxDict = taxonomy_dictionary("rankedlineage_tabRemoved.dmp")

    # if taxonomy id from blast file exist in the dictionary
    for i in filBlast:
        if i[2] in taxDict:

            f.write(str('\t'.join(taxDict[i[2]]) + '\t' + '\t'.join(i) + '\n'))

    return


# The link_OTUtable function will store header and values of an OTU table
def link_OTUtable(OTUtable):
    tableDict = {}
    lableDict = {}

    with open(OTUtable, "r") as file:
        for line in file:
            ll = line.strip().split('\t')
            if "#" in ll[0]:
                lableDict[ll[0]] = ll[1:]  # storing header in dictionary
            else:
                # storing abundance values in a dictionary
                tableDict[ll[0]] = ll[1:]

    return {"tableDict": tableDict, "lableDict": lableDict}
