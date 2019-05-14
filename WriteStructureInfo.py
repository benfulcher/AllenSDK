# ADAPTED FROM LayerStructuresIsocortex.ipynb
from allensdk.api.queries.ontologies_api import OntologiesApi
from allensdk.core.structure_tree import StructureTree
import pandas as pd
#-------------------------------------------------------------------------------

def writeIDsAndInfo(structureList,IDFileName,infoFileName):
    'Writes structure IDs and information to file'
    # Write IDs to .csv file:
    structIDs = [x['id'] for x in structureList]
    df = pd.DataFrame(structIDs)
    df.to_csv(IDFileName, index=False, header=False)

    # Write layer info:
    structInfo = [[ind+1,x['id'],x['acronym'],x['name']] for ind,x in enumerate(structureList)]
    df = pd.DataFrame(structInfo)
    df.to_csv(infoFileName, index=False, header=False)
def getFullStructureTree():
    oapi = OntologiesApi()
    structure_graph = oapi.get_structures_with_sets([1])  # 1 is the id of the adult mouse structure graph
    # structure_graph = oapi.get_structures(structure_graph_ids=1) # (Alternative method)
    # This removes some unused fields returned by the query
    structure_graph = StructureTree.clean_structures(structure_graph)
    tree = StructureTree(structure_graph)
    return tree
def writeAllCorticalStructures(tree):
    nameMap = tree.get_name_map()
    cortexStructures = [name_map[i] for i in tree.descendant_ids([315])[0]]
def writeLayerStructures(tree):
    'Write layer-related cortical structures'
    nameMap = tree.get_name_map()
    # Layer-related isocortical structures:
    layerIDs = [i for i in tree.descendant_ids([315])[0] if 'layer' in name_map[i] or 'Layer' in name_map[i]]
    layerStructs = tree.get_structures_by_id(layerIDs)
    writeIDsAndInfo(layerStructs,'structIDs_layers.csv','structInfo_layers.csv')
def writeCorticalStructures(tree):
    nameMap = tree.get_name_map()
    # Non-layer-related children of isocortex:
    notLayerIDs = [i for i in tree.descendant_ids([315])[0] if 'layer' not in name_map[i] and 'Layer' not in name_map[i] and name_map[i]!='Isocortex']
    notLayerStructs = tree.get_structures_by_id(notLayerIDs)
    writeIDsAndInfo(notLayerStructs,'structIDs_cortex.csv','structInfo_cortex.csv')
def writePalomeroStructures(tree):
    theStructIDs = [1002,332,385,985,512,599,672,1080,165,718,733]
    palomeroStructs = tree.get_structures_by_id(theStructIDs)
    writeIDsAndInfo(palomeroStructs,'structIDs_palomero.csv','structInfo_palomero.csv')
#-------------------------------------------------------------------------------

tree = getFullStructureTree()
writeWhat = 'palomero'
if writeWhat=='palomero':
    writePalomeroStructures(tree)
elif writeWhat=='cortical':
    writeCorticalStructures(tree)
elif writeWhat=='layer':
    writeLayerStructures(tree)
