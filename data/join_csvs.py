
import glob
import pandas as pd
from sys import argv

args = argv

path = args[1] # use your path

allFiles = glob.glob(path + "/*.csv")
frame = pd.DataFrame()
list_ = []
for file_ in allFiles:
    df = pd.read_csv(file_, header=0)
    list_.append(df)
frame = pd.concat(list_, ignore_index = True )
