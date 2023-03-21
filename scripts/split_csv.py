import pandas as pd
import numpy as np

df = pd.read_table(snakemake.input.csv, sep=',',dtype={'subjects':str})
#df = pd.read_csv(snakemake.input.csv)
tsvs = np.array_split(df, df.shape[0])

for i, tsv in enumerate(tsvs):
    tsv.to_csv(snakemake.output.tsv.format(i), sep='\t', index=False, header=True)

#how to make (i) the subject ID?
