import pandas as pd

# Load in downloaded gene info (`AllGenes.py`):
geneInfoFilename = 'geneInfo.csv'
df = pd.read_csv(geneInfoFilename)
df.head()
# Get list of entrezIDs

geneList = ['Gabra1','Gabra2','Gabra3','Gabra4','Gabra5','Gabrb1','Gabbr1',
                'Gabbr2','Gria1','Gria2','Gria3','Gria4','Grin1','Grin2a',
                'Grin2b','Grin2c','Grin2d','Grin3a','Grin3b','Grina','Htr1a',
                'Htr1b','Htr1d','Htr1f','Htr2a','Htr2b','Htr2c','Htr3a','Htr3b',
                'Htr4','Htr5a','Htr5b','Htr6','Htr7','Chrm1','Chrm2','Chrm3']

entrezIDs = [df[df.acronym==gene].entrez_id.values[0] for gene in geneList]
print('Matched %u/%u genes' % (len(entrezIDs),len(geneList)))
dfEntrez = pd.DataFrame(entrezIDs)
dfEntrez.to_csv('entrezIDs_palomero.csv',index=False,header=False)
