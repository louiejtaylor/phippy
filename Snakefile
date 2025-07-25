import subprocess
from collections import defaultdict as dd

#include: rules/qc.smk
#include: rules/map.smk
#include: rules/summary.smk

SAMPLES = [str(s) for s in config["samples"].keys()]
DATA_DIR =str(config["io"]["data"])
OUTPUT_DIR = str(config["io"]["output"])

rule symlink_data:
    output:
        r1 = str(OUTPUT_DIR+"/qc/raw/{sample}_1.fastq.gz"),
        r2 = str(OUTPUT_DIR+"/qc/raw/{sample}_2.fastq.gz")
    run:
        in_r1 = str(DATA_DIR+"/"+readmap[wildcards.sample][0])
        in_r2 = str(DATA_DIR+"/"+readmap[wildcards.sample][0])
        # print(in_r1, in_r2)
        subprocess.run(["ln", "-sr", in_r1, output.r1])
        subprocess.run(["ln", "-sr", in_r2, output.r2])

rule all_symlink:
    input: expand(str(OUTPUT_DIR+"/qc/raw/{sample}_1.fastq.gz"), sample=SAMPLES)

rule uncompress_fastqs:
    input:
        gzipped = str(OUTPUT_DIR+"/qc/raw/{sample}_1.fastq.gz")
    output:
        unzipped = str(OUTPUT_DIR+"/qc/raw/{sample}_1.fastq")
    shell:
        """
        gunzip {input.gzipped} -c > {output.unzipped}
        """

rule extract_first_50:
    input:
        r1 = str(OUTPUT_DIR+"/qc/raw/{sample}_1.fastq")
    output:
        s50 = str(OUTPUT_DIR+"/summary/counts/samples/{sample}_first50.csv")
    threads: 1
    run:
        from Bio import SeqIO
        scount = dd(int)
        for record in SeqIO.parse(input.r1, format = "fastq"):
            scount[record.seq[:50]] += 1
        o = open(output.s50,"w")
        for k in scount.keys():
            o.write(str(scount[k])+","+str(k)+"\n")
        o.close()

rule all_extract:
    input: 
        expand(str(OUTPUT_DIR+"/summary/counts/samples/{sample}_first50.csv"), sample=SAMPLES)

rule combine_extracted:
    input:
        counts = expand(str(OUTPUT_DIR+"/summary/counts/samples/{sample}_first50.csv"), sample=SAMPLES)
    output:
        summary = str(OUTPUT_DIR+"/summary/counts/all_first50.csv")
    run:
        import pandas
        def merge_all(df_list):
            if len(df_list) == 1:
                return(df_list[0])
            else:
                return(pandas.merge(df_list[0], merge_all(df_list[1:]), on='seq', how='outer'))
        all_dfs = []
        # TODO: improve memory use; ONLY read in a df when merging it (not have all in memory)
        for i in input.counts:
            sample_name = i.split("/")[-1].replace("_first50.csv", "")
            all_dfs.append(pandas.read_csv(i, names = [sample_name, "seq"]))
        summary_df = merge_all(all_dfs)
        summary_df.to_csv(output.summary, index=False)

# make more elegant by having a single combination rule for all types of analyses (e.g. exact match vs. aligned)

rule incorporate_metadata:
    input:
        summary = str(OUTPUT_DIR+"/summary/counts/all_first50.csv"),
        md = str(config['metadata'])
    output:
        md_summary = str(OUTPUT_DIR+"/summary/counts/all_md_first50.csv")
    run:
        import pandas
        # TODO: check if md df has a 'seq' column 
        md_df = pandas.merge(pandas.read_csv(input.summary), pandas.read_csv(input.md), on='seq', how='outer')
        md_df.to_csv(output.md_summary, index=False)
