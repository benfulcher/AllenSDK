import pandas as pd
from allensdk.api.queries.rma_api import RmaApi
import allensdk.core.json_utilities as json_utilities
import numpy as np
import csv # For saving string data to csv

#-------------------------------------------------------------------------------
## Set filenames to output to:
json_file_name = 'allSections.json'
sectionDatasetFilename = 'sectionDatasetInfo.csv'
geneInfoFilename = 'geneInfo.csv'
entrezIDFilename = 'geneEntrezID.csv'

#-------------------------------------------------------------------------------
def GetAllSections():
    criteria = "[failed$eqFalse],products[id$eq1]"

    rows = QueryAPI(model='SectionDataSet',criteriaString=criteria,
                    includeString='genes',optionsString='[only$eq''genes.entrez_id,data_sets.id'']')
    return rows
#-------------------------------------------------------------------------------
def QueryAPI(model,criteriaString,includeString="",optionsString="",writeOut=[]):
    # Initiate RMA API for Allen data retrieval
    api = RmaApi()
    # Settings for retrieval
    rows = []
    blockSize = 2000
    done = False
    startRow = 0
    # for i in range(0, total_rows, blockSize):

    while not done:
        print "Row %d, attempting to retrieve %d rows..." % (startRow, blockSize)

        tot_rows = len(rows)
        if len(includeString)==0:
            rows += api.model_query(model=model,
                                    criteria=criteriaString,
                                    options=optionsString,
                                    start_row=startRow,
                                    num_rows=blockSize)
        else:
            rows += api.model_query(model=model,
                                    criteria=criteriaString,
                                    include=includeString,
                                    options=optionsString,
                                    start_row=startRow,
                                    num_rows=blockSize)

        numRows = len(rows) - tot_rows # additional rows retrieved on running the query
        startRow += numRows

        print "%d rows retrieved." % numRows

        # Check if we're at the end of the road
        if numRows == 0 or numRows < blockSize:
            done = True

        # Write out the results as they come in, if requested
        if isinstance(writeOut, basestring):
            json_utilities.write(json_file_name, rows)
            print "Wrote to %s" % json_file_name

    return rows
#-------------------------------------------------------------------------------
def GetAllGenes():
    criteria = "[organism_id$eq2][reference_genome_id$eq486545752]"
    rows = QueryAPI(model='Gene',criteriaString=criteria)
    return rows
#-------------------------------------------------------------------------------
def to_dataframe_genes(genes):
    fdata = []
    for gene in genes:
        fdata.append({
            'acronym': gene['acronym'],
            'id': gene['id'],
            'entrez_id': gene['entrez_id'],
        })
    return pd.DataFrame.from_records(fdata)
#-------------------------------------------------------------------------------
def SectionsToList(sections):
    # Given a list of sections, returns a list of entrez ids contained in those
    # section datasets
    entrezList = []
    for section in sections:
        if (len(section['genes'])>0):
            if isinstance(section['genes'][0]['entrez_id'],int):
                entrezList.append(section['genes'][0]['entrez_id'])
    entrezList.sort()
    return entrezList
#-------------------------------------------------------------------------------
def SaveSectionsToCSV(sections):
    # Given a list of sections, saves info about genes as csv, for reading in
    # to Matlab

    # Best is to have: 1) a section id to entrez id, and 2) data about all genes

    # 1. Section id to gene entrez id mapping (with plane of section also included):
    sectionList = []
    for section in sections:
        if (len(section['genes'])>0) and isinstance(section['genes'][0]['entrez_id'],int):
            sectionList.append({
                'section_id': section['id'],
                'entrez_id': section['genes'][0]['entrez_id'],
                'plane_of_section_id': section['plane_of_section_id']
                })
    # To dataframe:
    df_sections = pd.DataFrame.from_records(sectionList) #,index='section_id')
    df_sections = df_sections.sort('section_id')

    # Save as a .csv file:
    df_sections.to_csv(sectionDatasetFilename)

    # 2. All info about genes in one place
    geneList = []
    for section in sections:
        if (len(section['genes'])>0) and isinstance(section['genes'][0]['entrez_id'],int):
            geneList.append({
                'acronym': section['genes'][0]['acronym'],
                'name': section['genes'][0]['name'],
                'gene_id': section['genes'][0]['id'],
                'entrez_id': section['genes'][0]['entrez_id']
                })
    # To dataframe:
    df_genes = pd.DataFrame.from_records(geneList) #,index='entrez_id')
    df_genes = df_genes.sort('entrez_id')
    numGenesFull = df_genes.shape[0]
    df_genes = df_genes.drop_duplicates()
    numGenesFiltered = df_genes.shape[0]
    print "Genes filtered from %u to %u" % (numGenesFull, numGenesFiltered)
    # Save as a .csv file:
    df_genes.to_csv(geneInfoFilename)
#-------------------------------------------------------------------------------
def SaveListCSV(stringList,fileName):
    # Outputs a csv from a given list of strings
    resultFile = open(fileName,'wb')
    wr = csv.writer(resultFile, dialect='excel')
    wr.writerow(stringList)

#-------------------------------------------------------------------------------

# --1-- Retrieve and save section information:
# Download and save all of the gene data to file:
sections = GetAllSections()
# Save section data to csv:
SaveSectionsToCSV(sections)
# Saves to:
# - sectionDatasetFilename
# - geneInfoFilename

# --2-- Get gene information:
# Get unique entrez IDs
geneEntrezList = SectionsToList(sections)
entrezSet = set(geneEntrezList)
geneEntrezList = list(entrezSet)
geneEntrezList.sort()
print "There are %u unique genes in section datasets" % len(entrezSet)
SaveListCSV(geneEntrezList,entrezIDFilename)

# Saves to:
# - entrezIDFilename


# genes = GetAllGenes()
# df = df.sort('entrez_id')
# df.entrez_id.unique
# df = to_dataframe(genes)
# df = df.sort('entrez_id')
# df.entrez_id.unique
# gb = df.groupby('structure_id')
