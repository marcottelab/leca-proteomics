import pandas as pd
import sys
import os
import time
from tqdm lsimport tqdm

start_time = time.time()
writer = pd.ExcelWriter('leca_ppis.xlsx') # Arbitrary output name
for csvfilename in tqdm(sys.argv[1:]):
    df = pd.read_csv(csvfilename)
    df.to_excel(writer,sheet_name=os.path.splitext(csvfilename)[0])
writer.save()
print("Runtime = {} seconds.".format(time.time() - start_time))