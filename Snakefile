import subprocess, pandas
from collections import defaultdict as dd

#include: rules/qc.smk
#include: rules/map.smk
#include: rules/summary.smk

from pathlib import Path

SAMPLE_MAP = config["samples"]
SAMPLES = [str(s) for s in SAMPLE_MAP.keys()]

DB_MAP = config["dbs"]
DBS = [str(s) for s in DB_MAP.keys()]

DATA_DIR = str(config["io"]["data"])
OUTPUT_DIR = str(config["io"]["output"])


summary_n = config["io"]["n_bases"]

rule symlink_data:
    output:
        r1 = str(OUTPUT_DIR+"/qc/raw/{sample}_1.fastq.gz"),
        r2 = str(OUTPUT_DIR+"/qc/raw/{sample}_2.fastq.gz")
    run:
        in_r1 = str(DATA_DIR+"/"+SAMPLE_MAP[wildcards.sample][0])
        in_r2 = str(DATA_DIR+"/"+SAMPLE_MAP[wildcards.sample][1])
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

rule symlink_peptide_map:
    output:
        pep_map = str(OUTPUT_DIR+"/dbs/pep_maps/{db}.csv")
    run:
        in_map = str(DB_MAP[wildcards.db]["fp"])
        subprocess.run(["ln", "-sr", in_map, output.pep_map])

rule map_to_peptides:
    input:
        pep_summary = str(OUTPUT_DIR+"/summary/counts/raw/{sample}_first"+str(summary_n)+".csv"),
        pep_db = str(OUTPUT_DIR+"/dbs/pep_maps/{db}.csv")
    output:
        mapped = str(OUTPUT_DIR+"/summary/counts/mapped/{db}/{sample}_first"+str(summary_n)+".csv")
    params:
        sample = "{sample}",
        db = "{db}",
        n = int(summary_n)
    run:
        from math import floor
        # keys: format id_col seq_col
        dbinfo = DB_MAP[params.db]
        trunc_n = params.n
        samp_join_on = "seq"
        if dbinfo["format"] == "aa":
            trunc_n = floor(trunc_n / 3)
            samp_join_on = "aas"
        map_df = pandas.read_csv(input.pep_db)
        map_df["_seq_n"] = map_df.apply(lambda r: r[dbinfo["seq_col"]][:trunc_n], axis = 1)
        mapped_df = pandas.merge(pandas.read_csv(input.pep_summary,names = [params.sample, "seq", "aas"]), 
                                 map_df[["_seq_n",dbinfo["id_col"]]], 
                                 left_on=samp_join_on, right_on='_seq_n', how='outer'
                                )
        collapsed_summary = mapped_df[[dbinfo["id_col"],params.sample]].groupby(dbinfo["id_col"],dropna=False).sum()
        collapsed_summary.to_csv(output.mapped)

# for large datasets, combine_extracted has a large memory requirement. options to mitigate:
# (1) map on a samplewise basis. low memory, but repeated computations (probably fine)
# (2) analyze only peptide sequence (not nt)
# (3) filter out peptides with internal stops

rule summarize_mapped:
    input:
        counts = expand(str(OUTPUT_DIR+"/summary/counts/mapped/{{db}}/{sample}_first"+str(summary_n)+".csv"), sample=SAMPLES)
    output:
        summary = str(OUTPUT_DIR+"/summary/counts/mapped/{db}_all_mapped_first"+str(summary_n)+".csv")
    params:
        n = int(summary_n),
        db = "{db}"
    run:
        def _open_merge_all(fp_list):
            """
            Helper function to recursively open and merge multiple counts tables.
            """
            dbinfo = DB_MAP[params.db]
            # sample_name = fp_list[0].split("/")[-1].replace("_first"+str(summary_n)+".csv", "")
            if len(fp_list) == 1:
                return(pandas.read_csv(fp_list[0]))
            else:
                new_df = pandas.read_csv(fp_list[0])
                return(pandas.merge(new_df, _open_merge_all(fp_list[1:]), on=dbinfo["id_col"], how='outer'))

        summary_df = _open_merge_all(input.counts)
        summary_df.to_csv(output.summary, index=False)

rule all_map:
    input:
        maps = expand(str(OUTPUT_DIR+"/summary/counts/mapped/{db}_all_mapped_first"+str(summary_n)+".csv"), db=DBS)

rule de_novo_annotate:
    input:
        summary = str(OUTPUT_DIR+"/summary/counts/raw/{sample}_first"+str(summary_n)+".csv"),
        db = str(config['io']['pep_fasta'])
    output:
        annotated = str(OUTPUT_DIR+"/summary/annotated/ann_{sample}_first"+str(summary_n)+".csv")
    threads: 1
    run:
        from Bio import SeqIO
        prot_db = [(str(record.seq), record.id) for record in SeqIO.parse(input.db, "fasta")]
        o = open(output.annotated, 'w')
        counter = 0
        with open(input.summary, 'r') as f:
            for line in f:
                pep = line.split(",")[2].strip()
                if "*" in pep[:-1]: # stop codon at last pos is fine
                    o.write(line.strip()+',internal_stop\n')
                    continue
                else:
                    if pep[-1] == "*":
                        pep = pep[:-1]

                matches = []
                for record in prot_db:
                    pos = record[0].find(pep)
                    if pos >= 0:
                        matches.append([pos,record[1]])
                if len(matches) == 0:
                    o.write(line.strip()+',no_exact_matches\n')
                else:
                    o.write(line.strip() + "," + ";".join(["~".join([str(i) for i in m]) for m in matches])+'\n')


        o.close()

rule build_map_denovo:
    input:
        annotated = expand(str(OUTPUT_DIR+"/summary/annotated/ann_{sample}_first"+str(summary_n)+".csv"), sample=SAMPLES)
    output:
        map = str(OUTPUT_DIR+"/summary/de_novo_map.csv")
    threads: 1
    run:
        # count,ntseq,aaseq,position~accession;pos2~acc2
        def _open_dedup(fp_list):
            """
            Helper function to recursively generate aa-accession map.
            Drops unannotated lines.
            """
            # correctly breaks if fp_list is empty
            df = pandas.read_csv(fp_list[0], header=None)
            df.drop(df[df[3] == "no_exact_matches"].index, inplace=True)
            df.drop(df[df[3] == "internal_stop"].index, inplace=True)
            df.drop(df.columns[[0,1]], axis=1, inplace=True)
            if len(fp_list) == 1:
                return(df.drop_duplicates())
            else:
                return(pandas.concat([df, _open_dedup(fp_list[1:])], ignore_index=True).drop_duplicates())

        _open_dedup(input.annotated).to_csv(output.map, index=False)

rule preprocess_annotation:
    input:
        annotated = str(OUTPUT_DIR+"/summary/annotated/ann_{sample}_first"+str(summary_n)+".csv")
    output:
        preprocessed = str(OUTPUT_DIR+"/summary/annotated/intermediates/summ_{sample}_first"+str(summary_n)+".csv")
    threads: 1
    params: sample = "{sample}"
    run:
        df = pandas.read_csv(input.annotated, header = None)
        df.drop(df.columns[[1,2]], axis=1, inplace=True)
        df.columns = [params.sample, "seqs"]
        collapsed_summary = df.groupby("seqs",dropna=False).sum()
        collapsed_summary.to_csv(output.preprocessed)


rule all_annotate:
    input:
        preprocessed = expand(str(OUTPUT_DIR+"/summary/annotated/intermediates/summ_{sample}_first"+str(summary_n)+".csv"), sample=SAMPLES)
    output:
        annotated_summary = str(OUTPUT_DIR+"/summary/annotated/all_first_"+str(summary_n)+"_hits_annotated.csv")
    threads: 1
    run:
        first = True
        for fp in input.preprocessed:
            df = pandas.read_csv(fp)
            if first:
                merged_df = df.copy()
                first = False
            else:
                merged_df = pandas.merge(df, merged_df.copy(), on="seqs", how='outer')
            del df
        merged_df.to_csv(output.annotated_summary, index=False)

