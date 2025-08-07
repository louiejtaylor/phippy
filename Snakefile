import subprocess
from collections import defaultdict as dd

#include: rules/qc.smk
#include: rules/map.smk
#include: rules/summary.smk

SAMPLES = [str(s) for s in config["samples"].keys()]
DATA_DIR =str(config["io"]["data"])
OUTPUT_DIR = str(config["io"]["output"])

SAMPLE_MAP = config["samples"]

rule symlink_data:
    output:
        r1 = str(OUTPUT_DIR+"/qc/raw/{sample}_1.fastq.gz"),
        r2 = str(OUTPUT_DIR+"/qc/raw/{sample}_2.fastq.gz")
    run:
        in_r1 = str(DATA_DIR+"/"+SAMPLE_MAP[wildcards.sample][0])
        in_r2 = str(DATA_DIR+"/"+SAMPLE_MAP[wildcards.sample][0])
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
        def open_merge_all(fp_list):
            """
            Helper function to recursively open and merge multiple counts tables.
            """
            sample_name = fp_list[0].split("/")[-1].replace("_first50.csv", "")
            if len(fp_list) == 1:
                return(pandas.read_csv(fp_list[0], names = [sample_name, "seq"]))
            else:
                new_df = pandas.read_csv(fp_list[0], names = [sample_name, "seq"])
                return(pandas.merge(new_df, open_merge_all(fp_list[1:]), on='seq', how='outer'))

        summary_df = open_merge_all(input.counts)
        summary_df.to_csv(output.summary, index=False)



# make more elegant by having a single combination rule for all types of analyses (e.g. exact match vs. aligned)

rule incorporate_metadata:
    input:
        summary = str(OUTPUT_DIR+"/summary/counts/all_first50.csv"),
        md = str(config['io']['metadata'])
    output:
        md_summary = str(OUTPUT_DIR+"/summary/counts/all_md_first50.csv")
    run:
        import pandas
        # TODO: check if md df has a 'seq' column 
        md_df = pandas.merge(pandas.read_csv(input.summary), pandas.read_csv(input.md), left_on='seq', right_on="oligo_first_50", how='outer')
        md_df.to_csv(output.md_summary, index=False)


rule summarize_collapse_unmapped:
    input:
        md_summary = str(OUTPUT_DIR+"/summary/counts/all_md_first50.csv")
    output:
        md_collapsed = str(OUTPUT_DIR+"/summary/counts/all_md_collapsed_first50.csv")
    run:
        import pandas
        df_summary = pandas.read_csv(input.md_summary)
        collapsed_summary = df_summary[["pep_id"]+SAMPLES].groupby("pep_id",dropna=False).sum()
        collapsed_summary.to_csv(output.md_collapsed)

rule extract_translate_first_120:
    input:
        r1 = str(OUTPUT_DIR+"/qc/raw/{sample}_1.fastq")
    output:
        s120 = str(OUTPUT_DIR+"/summary/counts/translations/{sample}_first120.csv")
    threads: 1
    run:
        from Bio import SeqIO
        from Bio.Seq import Seq
        scount = dd(int)
        for record in SeqIO.parse(input.r1, format = "fastq"):
            scount[record.seq[:120]] += 1
        o = open(output.s120,"w")
        for k in scount.keys():
            o.write(str(scount[k])+","+str(Seq(k).translate())+","+str(k)+"\n")
        o.close()

rule all_extract_120:
    input:
        expand(str(OUTPUT_DIR+"/summary/counts/translations/{sample}_first120.csv"), sample=SAMPLES)
