rule reference_download:
  input:
    fastqdir="resources/fastq/{sample_fbc}",
    libraries="config/samples/fbcs/{sample_fbc}_libraries.csv",
    features="config/samples/fbcs/{sample_fbc}_features.csv",