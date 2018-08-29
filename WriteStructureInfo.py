# ADAPTED FROM LayerStructuresIsocortex.ipynb
from allensdk.api.queries.ontologies_api import OntologiesApi
from allensdk.core.structure_tree import StructureTree
import pandas as pd
#-------------------------------------------------------------------------------

def writeIDsAndInfo(structureList,IDFileName,infoFileName):
    # Write IDs to .csv file:
    structIDs = [x['id'] for x in structureList]
    df = pd.DataFrame(structIDs)
    df.to_csv(IDFileName, index=False, header=False)

    # Write layer info:
    structInfo = [[ind+1,x['id'],x['acronym'],x['name']] for ind,x in enumerate(structureList)]
    df = pd.DataFrame(structInfo)
    df.to_csv(infoFileName, index=False, header=False)

#-------------------------------------------------------------------------------
oapi = OntologiesApi()
structure_graph = oapi.get_structures_with_sets([1])  # 1 is the id of the adult mouse structure graph
# structure_graph = oapi.get_structures(structure_graph_ids=1) # (Alternative method)
# This removes some unused fields returned by the query
structure_graph = StructureTree.clean_structures(structure_graph)
tree = StructureTree(structure_graph)

# Get names of all structures in the ABA Isocortex that contain the word 'layer':
name_map = tree.get_name_map()
# cortexStructures = [name_map[i] for i in tree.descendant_ids([315])[0]]
# layerStructures = [k for k in a if 'layer' in k]

#-------------------------------------------------------------------------------
# Layer-related isocortical structures:
layerIDs = [i for i in tree.descendant_ids([315])[0] if 'layer' in name_map[i] or 'Layer' in name_map[i]]
layerStructs = tree.get_structures_by_id(layerIDs)
writeIDsAndInfo(layerStructs,'structIDs_layers.csv','structInfo_layers.csv')

#-------------------------------------------------------------------------------
# Non-layer-related children of isocortex:
notLayerIDs = [i for i in tree.descendant_ids([315])[0] if 'layer' not in name_map[i] and 'Layer' not in name_map[i] and name_map[i]!='Isocortex']
notLayerStructs = tree.get_structures_by_id(notLayerIDs)
writeIDsAndInfo(notLayerStructs,'structIDs_cortex.csv','structInfo_cortex.csv')


# ABAIsocortex = [tree.get_structures_by_id([i])[0]['acronym'] for i in tree.descendant_ids([315])[0]]
# df = pd.DataFrame(ABAIsocortex)
# df.to_csv('ABAIsocortex.csv', index=False, header=False)
