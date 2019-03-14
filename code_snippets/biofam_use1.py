from biofam.run.entry_point import entry_point

infiles = ["gene_exp_case.txt", "methylation_case.txt",
           "gene_exp_control.txt", "methylation_control.txt"]

feature_groups =  ["expression", "methylation","expression", "methylation"]
sample_groups = ["case", "case", "control", "control"]
lik = ["gaussian", "bionomial"]

outfile = dir+"results.hdf5"

ent = entry_point()
ent.set_data_options(lik)
ent.set_data_from_files(infiles, views, groups)
ent.set_model_options(ard_z=True, sl_w=False, sl_z=False, ard_w=True, noise_on='features')
ent.set_train_options()
ent.build()
ent.run()
ent.save(outfile)
