import pandas as pd
from allensdk.api.queries.rma_api import RmaApi
import allensdk.core.json_utilities as json_utilities
import numpy as np
import time
import csv # For saving string data to csv

#-------------------------------------------------------------------------------
# Globals:
#-------------------------------------------------------------------------------
# graph_id = 1 is the "Mouse Brain Atlas"
# (cf. http://help.brain-map.org/display/api/Atlas+Drawings+and+Ontologies)
mouse_graph_id = 1
mouse_product_id = 1
json_file_name = 'unionizes.json'

def download_mouse_unionizes(geneEntrezID,structureIDs,doReduced=False,writeOut=False):
    # Given genes and set of structures, retrieves expression data from Allen API
    # Can set geneEntrezID to empty to retrieve for all genes
    # Can set geneEntrezID to be a single number (array length 1) for a single gene
    # Can set geneEntrezID to an array to retrieve for multiple genes
    # doReduced aims to speed up the process by only retrieving bare-bones info

    # Convert list of structure IDs to a comma-separated string:
    string_structureIDs = ','.join(str(sid) for sid in structureIDs)

    # Set defaults:
    if doReduced:
        print("Only retrieving a reduced subset of results from the queries")
        optionsString = "[only$eq'structures.id,data_sets.id,data_sets.plane_of_section_id,genes.entrez_id,structure_unionizes.expression_energy,structure_unionizes.expression_density']"
    else:
        optionsString = None

    # Set criteria :
    if len(geneEntrezID)==0:
        # Get for all genes
        print("Get for all genes")
        # cf. Full query: http://api.brain-map.org/api/v2/data/query.xml?
        # criteria=model::StructureUnionize,rma::criteria,structure[id$eq22][graph_id$eq1],
        # section_data_set[failed$eqFalse](products[id$eq1]),
        # rma::include,section_data_set,structure,section_data_set(genes)
        criteria = [
            "structure[id$in%s][graph_id$eq%d]," % (string_structureIDs, mouse_graph_id),
            "section_data_set[failed$eqfalse](products[id$eq%d])" % (mouse_product_id)
            ]
    elif len([geneEntrezID])==1:
        # Get for a single gene:
        print("Get for the single gene: entrez_id = %u" % geneEntrezID[0])
        criteria = [
            "structure[id$in%s][graph_id$eq%d]," % (string_structureIDs, mouse_graph_id),
            "section_data_set[failed$eqfalse](products[id$eq%d],genes[entrez_id$eq%d])" % (mouse_product_id, geneEntrezID[0])
            ]
    elif len(geneEntrezID)>1:
        # Get for multiple genes:
        print("Getting expression results for %u genes" % len(geneEntrezID))
        string_geneEntrezIDs = ','.join(str(eid) for eid in geneEntrezIDs)
        criteria = [
            "structure[id$in%s][graph_id$eq%d]," % (string_structureIDs, mouse_graph_id),
            "section_data_set[failed$eqfalse](products[id$eq%d],genes[entrez_id$in%s])" % (mouse_product_id, string_geneEntrezIDs)
            ]

    # Query the API:
    expressionData = QueryAPI(model='StructureUnionize',criteriaString="".join(criteria),
                includeString="section_data_set,structure,section_data_set(genes)",
                optionsString=optionsString,writeOut=writeOut)

    return expressionData

def QueryAPI(model,criteriaString,includeString=None,optionsString=None,writeOut=False):
    # Send a query to the Allen API, and assemble results

    # Initiate RMA API for Allen data retrieval
    api = RmaApi()

    # Settings for retrieval
    rows = []
    blockSize = 2000
    done = False
    startRow = 0
    # for i in range(0, total_rows, blockSize):

    while not done:
        print("Row %d, attempting to retrieve %d rows..." % (startRow, blockSize))

        tot_rows = len(rows)

        # apiQueryPartial = partial(api.model_query,model=model,criteria=criteriaString,
                                    # startRow=startRow,num_rows=blockSize)

        rows += api.model_query(model=model,
                                criteria=criteriaString,
                                include=includeString,
                                options=optionsString,
                                start_row=startRow,
                                num_rows=blockSize)

        numRows = len(rows) - tot_rows # additional rows retrieved on running the query
        startRow += numRows

        print("%d rows retrieved." % numRows)

        # Check if we're at the end of the road
        if numRows == 0 or numRows < blockSize:
            done = True

        # Write out the results as they come in, if requested
        if writeOut:
            json_utilities.write(json_file_name, rows)
            print("Wrote to %s" % json_file_name)

    return rows

def unionizes_to_dataframe(unionizes):
    fdata = []
    for unionize in unionizes:
        fdata.append({
            'structure_id': unionize['structure_id'],
            'expression_energy': unionize['expression_energy'],
            'expression_density': unionize['expression_density'],
            'data_set_id': unionize['section_data_set']['id'],
            'plane_of_section_id': unionize['section_data_set']['plane_of_section_id'],
            # 'gene_acronym': unionize['section_data_set']['genes'][0]['acronym'],
            # 'gene_name': unionize['section_data_set']['genes'][0]['name'],
            # 'gene_entrez': unionize['section_data_set']['genes'][0]['entrez_id'],
        })
    return pd.DataFrame.from_records(fdata)

def SaveListCSV(stringList,fileName):
    resultFile = open(fileName,'wb')
    wr = csv.writer(resultFile, dialect='excel')
    wr.writerow(stringList)

# def GetGeneDetails(entrezIDs,filterFields,graph_id=1):
#
#     [organism_id$eq2][reference_genome_id$eq486545752]

def GetStructureDetails(structureIDs,filterFields,graph_id=1):
    # Get structure info for structureIDs given from Allen API
    # (mouse brain atlas)

    print("\n---Downloading all mouse brain structures using the Allen Institute API...\n")

    # Convert list of structure IDs to a comma-separated string:
    string_structureIDs = ','.join(str(sid) for sid in structureIDs)

    # Get the structure info by querying Allen API:
    structs = QueryAPI(model="Structure",criteriaString=("[id$in%s][graph_id$eq%u]" \
                                % (string_structureIDs,graph_id))) # (Mouse Brain Atlas)

    print("%u mouse-brain structures retrieved!\n" % len(structs))

    # Filter to only certain fields, specified as filterFields:
    structsFilt = []
    for i,s in enumerate(structs):
        structsFilt.append({keep_key: s[keep_key] for keep_key in filterFields})

    return structsFilt

    # # Build a dict from structure id to structure and identify each node's direct descendants:
    # structHash = {} # dict
    # for s in structs:
    #     s['num_children'] = 0 # add this field
    #     s['structure_id_path'] = [int(sid) for sid in s['structure_id_path'].split('/') if sid != '']
    #     s['acronym'] = s['acronym'].rstrip() # some acronyms have trailing whitespace (e.g., 'SUM'); remove
    #     structHash[s['id']] = s
    #
    # # Iterate through to add the 'num_children' field by incrementing the 'num_children'
    # # field given structures with parents:
    # for sid,s in structHash.iteritems():
    #     if len(s['structure_id_path']) > 1:
    #         parentId = s['structure_id_path'][-2]
    #         structHash[parentId]['num_children'] += 1


def main():
    # Download and save all of the data to file:

# INPUT file names:
structIDSource = 'ValerioStructIDs.csv' # 'layerIDs.csv' # 'structureIDs.csv'
entrezSource = 'ValerioEntrezIDs.csv' # 'geneEntrezID.csv'

# OUTPUT file names:
structInfoFilename = 'structureInfoAdra.csv' # structureInfo.csv
allDataFilename = 'bigDataFrameAdra.csv'

#---------------------------------------------------------------------------
# Read in structure IDs from csv -> save detailed structure info to file
#---------------------------------------------------------------------------
structureIDs = np.genfromtxt(structIDSource,delimiter=',')
print("Retrieved %d structures from %s..." % (len(structureIDs),structIDSource))
# Get details of the structures:
filterFields = ['id','name','acronym','color_hex_triplet']
structs = GetStructureDetails(structureIDs,filterFields,1)
# To dataframe:
df_struct = pd.DataFrame.from_records(structs)
# Save as a .csv file:
df_struct.to_csv(structInfoFilename)

#---------------------------------------------------------------------------
# RETRIEVING GENES ONE AT A TIME
#---------------------------------------------------------------------------
# Kind of need to go one at a time since the request is too long if you try
# to get all genes at once?

# Read in gene entrez IDs from csv:
geneEntrezIDs = np.genfromtxt(entrezSource,delimiter=',')
print("Read in %d genes from %s..." % (len(geneEntrezIDs), entrezSource))

# geneEntrezIDs = np.concatenate((geneEntrezIDs[0:1],[20604.]),axis=0)
unionizes = []
startTime = time.time()
for ind, geneEntrezID in enumerate(geneEntrezIDs):
    unionizes_gid = download_mouse_unionizes([geneEntrezID],structureIDs,
                                            doReduced=True,writeOut=False)
    unionizes = unionizes + unionizes_gid
    # df_i = unionizes_to_dataframe(unionizes_gid)
    # df = pd.concat([df,df_i])
    timeSoFar = time.time() - startTime
    timeRemaining = timeSoFar/(ind+1)*(len(geneEntrezIDs)-ind+1)/60
    print("We're at %d/%d, %f min remaining" % (ind+1,len(geneEntrezIDs),timeRemaining))
totalTime = (time.time() - startTime)/60.
print("Took %f min total" % totalTime)

# Make a data frame of the datasets downloaded:
df = unionizes_to_dataframe(unionizes)
df = df.sort_values(by=['structure_id','data_set_id'])
df.to_csv(allDataFilename)
print("Saved data frame to %s" % allDataFilename)
print(df)

# df = pd.read_csv(allDataFilename)

#-------------------------------------------------------------------------------
# RETRIEVING GENES ALL TOGETHER IN ONE QUERY
#-------------------------------------------------------------------------------
# Another option would be to retrieve everything in one query, and then filter
# back to the matching genes...? (25924 rows per structure...):
# No need to write to a json file as you go (slows down)
# Only retrieve the necessary fields (doReduced)
# unionizes = download_mouse_unionizes([],structureIDs,doReduced=True,writeOut=False)

#-------------------------------------------------------------------------------
# Write the agglomerated datasets to file
# json_utilities.write(json_file_name, unionizes)
# print("Wrote %d unionizes to %s" % (len(unionizes),json_file_name))

# Read in previously-downloaded results:
# unionizes = json_utilities.read(json_file_name)

# To dataframe and then filter to include only genes in our geneEntrezIDs list:
# dfFull = unionizes_to_dataframe(unionizes)
# dfFull.shape[0]
# # Now we need to filter to include only genes in our geneEntrezIDs list:
# # (actually this filtering is not necessary) -- potentially a nice check that we only retrieved
# # entrez_ids corresponding to section datasets with expression (22,000)
# df = dfFull.copy()
# df = df[df['gene_entrez'].isin(geneEntrezIDs)]
# df.shape[0]
# len(geneEntrezIDs)
# # (so something like 3,000 duplicates?)
# df = df.sort(['structure_id','data_set_id'])
# df.tail
# df_filter = df.drop_duplicates()
# df['data_set_id']

# # Save out just the gene information:
# dfGenes = df[]
#
# # Save data summary as a .csv file:
# csv_file_name = 'filteredDataFrame.csv'
# df.to_csv(csv_file_name)

# Compute the mean expression energy for each structure_id:
# gb = df.groupby(['gene_entrez','structure_id'])['expression_energy'].mean()
# print(gb)

# gb.pivot(index="structure_id", columns="gene_entrez", values="expression_energy")
print("Writing expression energy and density matrices out")
for expressionType in ['expression_energy','expression_density']:
    table = pd.pivot_table(df, values=expressionType, index="structure_id",
                        columns="data_set_id", aggfunc=np.mean)
    # table = pd.pivot_table(df, values='expression_energy', index="structure_id",
    #                     columns="gene_entrez", aggfunc=np.mean)
    print(table)

    # Save mean expression energy to a .csv output file:
    fileName = ("%s_%ux%u.csv" % (expressionType, table.shape[0], table.shape[1]))
    table.to_csv(fileName,na_rep='NaN')
    print("Saved expression energy of %u datasets over %u structures to %s" \
                    % (table.shape[1],table.shape[0],fileName))

# Need to also output the section dataset IDs (columns) and struct IDs (rows)
dataSetIDs = pd.DataFrame(list(table.axes[0]))
dataSetIDs.to_csv('structIDs_Rows.csv', index=False, header=False)
dataSetIDs = pd.DataFrame(list(table.axes[1]))
dataSetIDs.to_csv('dataSetIDs_Columns.csv', index=False, header=False)

    # Summary to screen
    # gdf = gb.agg({'data_set_id': pd.Series.nunique})
    # print(gdf)

if __name__ == "__main__": main()
