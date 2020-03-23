
from mofapy2.run.entry_point import entry_point
from re import search
import pandas as pd
import numpy as np
import os

################
## Define I/O ##
################

host = os.uname()[1]
if search("ricard", host):
    datadir = '/Users/ricard/data/mofa2_simulations/test'
    outfile = '/Users/ricard/data/mofa2_simulations/test/hdf5/test.hdf5'
elif search("embl", host):
	datadir = "/g/stegle/ricard/mofa2_simulations/test"
	outfile = "/g/stegle/ricard/mofa2_simulations/test/hdf5/test.hdf5"
else:
    print("Computer not recognised"); exit()

###############
## Load data ##
###############

datafile = "%s/data.txt.gz" % datadir
data = pd.read_csv(datafile, delimiter="\t", header=0)

##########
## MOFA ##
##########

# initialise
ent = entry_point()

# Set data (data.frame format)
ent.set_data_df(data)

# Set model options
ent.set_model_options(factors=15, spikeslab_factors=False, spikeslab_weights=False)

# Set training options
ent.set_train_options(iter=100, seed=42, convergence_mode="fast")

# Set stochastic inference options
# ent.set_stochastic_options(learning_rate=0.5, forgetting_rate=0.5, batch_size=0.5, start_stochastic=1)

# Build the model
ent.build()

# Train the model
ent.run()

# Save the model
ent.save(outfile)