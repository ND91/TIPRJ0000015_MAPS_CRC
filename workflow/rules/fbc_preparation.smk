rule fbc_libraries:
  input:
    config['sample_files']
  output:
    run_id_libraries="config/samples/fbcs/{run_fbc}_libraries.csv",
  message:
    "--- Preparing the Cellranger library files for {params.run_fbc} ---"
  threads: 
    1
  resources:
    mem_mb=2048,
  params: 
    run_fbc="{run_fbc}",
  run:
    import os
  
    run_id_df = sample_files_fbc[sample_files_fbc["Run_ID"] == params.run_fbc]
    #run_id_df['Well10X_str'] = run_id_df['Well10X'].astype(str)
    
    #sample_fbc = run_id_df[['Run_ID', 'Well10X_str']].agg('-'.join, axis=1).astype(str) + "-ADT"
    sample_fbc = run_id_df['Run_ID'].astype(str) + "-ADT"
    fastqs_fbc = os.getcwd() + '/resources/fastq/' + sample_fbc.astype(str)
    
    libraries_fbc = pd.DataFrame(list(zip(fastqs_fbc,
                                          sample_fbc, 
                                          run_id_df["FBC_type"])),
                                 columns = ['fastqs', 'sample', 'library_type'])
    
    #sample_gex = run_id_df[['Run_ID', 'Well10X_str']].agg('-'.join, axis=1)
    sample_gex = run_id_df['Run_ID']
    fastqs_gex= os.getcwd() + '/resources/fastq/' + sample_gex.astype(str)
    
    libraries_gex = pd.DataFrame(list(zip(fastqs_gex, 
                                          sample_gex)),
                                 columns = ['fastqs', 'sample'])
    libraries_gex['library_type'] = pd.Series(["Gene Expression" for x in range(len(libraries_gex.index))])
    
    libraries = pd.concat([libraries_fbc, libraries_gex], ignore_index=True)
    libraries = libraries.drop_duplicates()
    
    libraries.to_csv(output.run_id_libraries, index = False)
