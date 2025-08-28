import subprocess, pandas
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
        s_n = str(OUTPUT_DIR+"/summary/counts/raw/{sample}_first"+str(summary_n)+".csv")
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
        expand(str(OUTPUT_DIR+"/summary/counts/raw/{sample}_first+"+str(summary_n)+".csv"), sample=SAMPLES)

rule map_to_peptide_db:
    input:
        pep_summary = str(OUTPUT_DIR+"/summary/counts/raw/{sample}_first"+str(summary_n)+".csv"),
        pep_db = str(config['io']['metadata'])
    output:
        mapped = str(OUTPUT_DIR+"/summary/counts/mapped/{sample}_first"+str(summary_n)+".csv")
    params:
        sample = "{sample}",
        n = int(summary_n)
    run:
        map_df = pandas.read_csv(input.pep_db)
        map_df["_seq_n"] = map_df.apply(lambda r: r.seq[:params.n], axis = 1)
        mapped_df = pandas.merge(pandas.read_csv(input.pep_summary,names = [params.sample, "seq", "aas"]), 
                                 map_df[["_seq_n","uid"]], 
                                 left_on='seq', right_on='_seq_n', how='outer'
                                )
        collapsed_summary = mapped_df[["uid",params.sample]].groupby("uid",dropna=False).sum()
        collapsed_summary.to_csv(output.mapped)

# for large datasets, combine_extracted has a large memory requirement. options to mitigate:
# (1) map on a samplewise basis. low memory, but repeated computations (probably fine)
# (2) analyze only peptide sequence (not nt)
# (3) filter out peptides with internal stops

rule summarize_mapped:
    input:
        counts = expand(str(OUTPUT_DIR+"/summary/counts/mapped/{sample}_first"+str(summary_n)+".csv"), sample=SAMPLES)
    output:
        summary = str(OUTPUT_DIR+"/summary/counts/all_mapped_first"+str(summary_n)+".csv")
    params: n = int(summary_n)
    run:
        def open_merge_all(fp_list):
            """
            Helper function to recursively open and merge multiple counts tables.
            """
            sample_name = fp_list[0].split("/")[-1].replace("_first"+str(summary_n)+".csv", "")
            if len(fp_list) == 1:
                return(pandas.read_csv(fp_list[0]))
            else:
                new_df = pandas.read_csv(fp_list[0])
                return(pandas.merge(new_df, open_merge_all(fp_list[1:]), on='uid', how='outer'))

        summary_df = open_merge_all(input.counts)
        summary_df.to_csv(output.summary, index=False)

#rule map_to_protein_db:
#    input:
#        pep_tables = str(OUTPUT_DIR+"/summary/counts/all_first"+str(summary_n)+".csv"),
#        db = str(config['io']['pep_fasta'])
#    output:
#        annotated = str(OUTPUT_DIR+"/summary/counts/annotations/hits_annotated.csv")
#    threads: 1
#    run:
#        from Bio import SeqIO
#        prot_db = list(SeqIO.parse(input.db,"fasta"))
#        # danger: reads full file into memory, need a lot of memory for big files
#        o = open(output.annotated, 'w')
#        counter = 0
#        with open(input.pep_tables, 'r') as f:
#            for line in f:
#                pep = line.split(",")[1]
#                if "*" in pep:
#                    o.write(line.strip()+',internal_stop\n')
#                else:
#                    matches = []
#                    for record in prot_db:
#                        pos = record.seq.find(pep)
#                        if pos >= 0:
#                            matches.append([pos,record.id])
#                    if len(matches) == 0:
#                        o.write(line.strip()+',no_exact_matches\n')
#                    else:
#                        o.write(line.strip()+ ";".join(["~".join([str(i) for i in m]) for m in matches])+'\n')
#                counter += 1
#                if counter % 100 == 0:
#                    print("processed "+str(counter))
#        o.close()

