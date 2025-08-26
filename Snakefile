import subprocess
from collections import defaultdict as dd

#include: rules/qc.smk
#include: rules/map.smk
#include: rules/summary.smk

from pathlib import Path

SAMPLES = [str(s) for s in config["samples"].keys()]
DATA_DIR = str(config["io"]["data"])
OUTPUT_DIR = str(config["io"]["output"])

SAMPLE_MAP = config["samples"]

summary_n = config["io"]["n_bases"]

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

rule extract_translate_first_n:
    input:
        r1 = str(OUTPUT_DIR+"/qc/raw/{sample}_1.fastq")
    output:
        s_n = str(OUTPUT_DIR+"/summary/counts/samples/{sample}_first"+str(summary_n)+".csv")
    params: n = int(summary_n)
    threads: 1
    run:
        from Bio import SeqIO
        from Bio.Seq import Seq
        scount = dd(int)
        for record in SeqIO.parse(input.r1, format = "fastq"):
            scount[record.seq[:params.n]] += 1
        o = open(output.s_n,"w")
        for k in scount.keys():
            o.write(str(scount[k])+","+str(k)+","+str(Seq(k).translate())+"\n")
        o.close()

rule all_extract:
    input:
        expand(str(OUTPUT_DIR+"/summary/counts/samples/{sample}_first+"+str(summary_n)+".csv"), sample=SAMPLES)

rule combine_extracted:
    input:
        counts = expand(str(OUTPUT_DIR+"/summary/counts/samples/{sample}_first"+str(summary_n)+".csv"), sample=SAMPLES)
    output:
        summary = str(OUTPUT_DIR+"/summary/counts/all_first"+str(summary_n)+".csv")
    params: n = int(summary_n)
    run:
        import pandas
        def open_merge_all(fp_list):
            """
            Helper function to recursively open and merge multiple counts tables.
            """
            sample_name = fp_list[0].split("/")[-1].replace("_first"+str(summary_n)+".csv", "")
            if len(fp_list) == 1:
                return(pandas.read_csv(fp_list[0], names = [sample_name, "seq", "aas"]))
            else:
                new_df = pandas.read_csv(fp_list[0], names = [sample_name, "seq", "aas"])
                return(pandas.merge(new_df, open_merge_all(fp_list[1:]), on='seq', how='outer'))

        summary_df = open_merge_all(input.counts)
        summary_df.to_csv(output.summary, index=False)

# make more elegant by having a single combination rule for all types of analyses (e.g. exact match vs. aligned)

rule incorporate_metadata:
    input:
        summary = str(OUTPUT_DIR+"/summary/counts/all_first"+str(summary_n)+".csv"),
        md = str(config['io']['metadata'])
    output:
        md_summary = str(OUTPUT_DIR+"/summary/counts/all_md_first"+str(summary_n)+".csv")
    run:
        import pandas
        # TODO: check if md df has a 'seq' column 
        md_df = pandas.merge(pandas.read_csv(input.summary), pandas.read_csv(input.md), on='seq', how='outer')
        md_df.to_csv(output.md_summary, index=False)


rule summarize_collapse_unmapped:
    input:
        md_summary = str(OUTPUT_DIR+"/summary/counts/all_md_first"+str(summary_n)+".csv")
    output:
        md_collapsed = str(OUTPUT_DIR+"/summary/counts/all_md_collapsed_first"+str(summary_n)+".csv")
    run:
        import pandas
        df_summary = pandas.read_csv(input.md_summary)
        collapsed_summary = df_summary[["Barcode ID"]+SAMPLES].groupby("Barcode ID",dropna=False).sum()
        collapsed_summary.to_csv(output.md_collapsed)

rule map_to_protein_db:
    input:
        pep_tables = str(OUTPUT_DIR+"/summary/counts/all_first"+str(summary_n)+".csv"),
        db = str(config['io']['pep_fasta'])
    output:
        annotated = str(OUTPUT_DIR+"/summary/counts/annotations/hits_annotated.csv")
    threads: 1
    run:
        from Bio import SeqIO
        prot_db = list(SeqIO.parse(input.db,"fasta"))
        # danger: reads full file into memory, need a lot of memory for big files
        o = open(output.annotated, 'w')
        counter = 0
        with open(input.pep_tables, 'r') as f:
            for line in f:
                pep = line.split(",")[1]
                if "*" in pep:
                    o.write(line.strip()+',internal_stop\n')
                else:
                    matches = []
                    for record in prot_db:
                        pos = record.seq.find(pep)
                        if pos >= 0:
                            matches.append([pos,record.id])
                    if len(matches) == 0:
                        o.write(line.strip()+',no_exact_matches\n')
                    else:
                        o.write(line.strip()+ ";".join(["~".join([str(i) for i in m]) for m in matches])+'\n')
                counter += 1
                if counter % 100 == 0:
                    print("processed "+str(counter))
        o.close()

