#### 1000G+ TE annotation ####

# Config file
configfile: "config.yaml"

import csv
from pathlib import Path

# Read the MANIFEST file to map genome names/IDs to their fasta paths
config["samples"] = {}
with open(config["MANIFEST"]) as fh:
    for line in csv.DictReader(fh, delimiter="\t"):  # Fixed to use `delimiter="\t"` for TSV
        config["samples"][line["sample"]] = line["genome_fa"]

# Helper function to set the working directory prefix for paths
def W(path):
    """This sets the prefix for all input/output paths."""
    return f"{config['WORK_DIR']}/{path}"

# Debugging: Print the genomes dictionary for verification
print(config["samples"])

### GENERAL RULE ###
rule all:
    input:
        # sglibs=expand(W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.TElib.fa"), sample=config["samples"].keys())
        # sglibs_updt=expand(W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.TElib.fa.LTRupdated.wgn.fa"), sample=config["samples"].keys())
        # libraries=expand(W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.TElib.mincopyfiltered.fa"), sample=config["samples"].keys())
        annotations=expand(W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.TEanno.sgannot.gff3"), sample=config["samples"].keys())
        # min3copy=expand(W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.TElib.mincopyfiltered.fa"), sample=config["samples"].keys())
        # allTEslib=W("make_te_lib/allTEs.fa")
        # panTElib=W("make_te_lib/panTElib.vsearch.centroids.fa")
        # panTE_reannot=expand(W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.TEanno.panTEannot.gff3"), sample=config["samples"].keys())
        # postprocessed_annot=expand(W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.TEanno.panTEannot.gff3.cleaned.gff"), sample=config["samples"].keys())
####################


### SINGLE GENOME PIPELINE ###
rule generate_EDTA_single_genome_TElib:
    input:
        fasta=lambda wc: config["samples"][wc.sample]
    output:
        TElib=W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.TElib.fa"),
        intactlib=W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.intact.fa")
    log: 
        "logs/{sample}_generate_EDTA_single_genome_TElib.log"
    threads: config["resources"]["threads"]
    params:
        directory=W("single_genome_EDTA/{sample}")
    conda: "panTEome"
    shell:
        """
        scripts_for_snakepipe/panTEome_generate_sg_libs.sh {input} {threads} {params.directory} > {log} 2>&1
        """

rule add_fullsize_LTRRTs:
    input:
        TElib=W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.TElib.fa"),
        intactlib=W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.intact.fa")
    output:
        W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.TElib.LTRupdated.fa.wgn.fa")
    log: 
        "logs/{sample}_add_fullsize_LTRRTs.log"
    threads: config["resources"]["threads"]
    params:
        directory=W("single_genome_EDTA/{sample}"),
        cleanup=config['scripts']["cleanup_nested"]
    conda: "panTEome"
    shell:
        """
        python3 scripts_for_snakepipe/add_fullsizeLTRRTs.py {input.TElib} {input.intactlib} {params.directory} {params.cleanup} {threads} > {log} 2>&1
        """

rule EDTA_annotate:
    input:
        tes_updated=W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.TElib.LTRupdated.fa.wgn.fa"),
        genome=lambda wc: config["samples"][wc.sample]
    output:
        W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.TEanno.sgannot.gff3"),
        W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.anno/{sample}.scaffolds_contigs.fa.mod.EDTA.RM.out")
    log:
        "logs/{sample}_EDTA_annotate.log"
    threads: config["resources"]["threads"]
    params:
        directory=W("single_genome_EDTA/{sample}")
    conda: "panTEome"
    shell:
        """
        scripts_for_snakepipe/panTEome_reannotate_from_sglib.sh {params.directory} {input.genome} {threads} > {log} 2>&1
        """

### panTE PIPELINE ###

rule extract_min_n_fl_copies:
    input: 
        gff=W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.TEanno.sgannot.gff3"),
        RMout=W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.anno/{sample}.scaffolds_contigs.fa.mod.EDTA.RM.out")
    output:
        W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.TElib.fa.keep.list"),
        W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.TElib.mincopyfiltered.fa")
    log:
        "logs/{sample}_extr_minnflcop.log"
    params: 
        min_fl_copy=config["params"]["min_fl_copy"],
        find_flTE=config['scripts']["find_flTE"]
    conda: "panTEome"
    shell:
        """
        bash scripts_for_snakepipe/get_n_fl_copies.sh {input.RMout} {params.min_fl_copy} {params.find_flTE} > {log} 2>&1 
        """

rule concat_all_telibs:
    input:
        expand(W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.TElib.mincopyfiltered.fa"), sample=config["samples"].keys())
    output:
        W("make_te_lib/allTEs.fa")
    shell:
        "cat {input} > {output} && sleep 1"


rule make_panTElib:
    input:
        allTElib=W("make_te_lib/allTEs.fa"),
        curatedlib=config["extra_datasets"]["curatedlib"]
    output:
        W("make_te_lib/panTElib.vsearch.centroids.fa")
    log:
        "logs/make_panTElib.log"
    threads: config["resources"]["threads"]
    conda: "panTEome"
    shell:
        """
        export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so.2
        python3 scripts_for_snakepipe/make_panTElib.py --all_TElib {input.allTElib} --curatedlib {input.curatedlib} --threads {threads} > {log} 2>&1
        echo 'done'  >> {log} 2>&1
        """


### PANGENOME ANNOTATION PIPELINE ###

rule reannotate_from_panTElib:
    input:
        gff=W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.TEanno.sgannot.gff3"),
        genome=lambda wc: config["samples"][wc.sample],
        panTElib=W("make_te_lib/panTElib.vsearch.centroids.fa")
    output:
        W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.TEanno.panTEannot.gff3"),
        W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.anno/{sample}.scaffolds_contigs.fa.mod.EDTA.TEanno.bed")
    log:
        "logs/{sample}_EDTA_reannotate_from_panTElib.log"
    threads: config["resources"]["threads"]
    params:
        directory=W("single_genome_EDTA/{sample}")
    conda: "panTEome"
    shell:
        """
        scripts_for_snakepipe/panTEome_reannotate_from_panTElib.sh {input.genome} {input.panTElib} {params.directory} {threads} > {log} 2>&1 
        """

rule postprocess_panTE_annotation:
    input:
        gff=W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.TEanno.panTEannot.gff3"),
        bed=W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.anno/{sample}.scaffolds_contigs.fa.mod.EDTA.TEanno.bed"),
        trashgff=lambda wildcards: f"{config['extra_datasets']['trash_annot_dir']}/{wildcards.sample}.TRASH.45S.gff"
    output:
        W("single_genome_EDTA/{sample}/{sample}.scaffolds_contigs.fa.mod.EDTA.TEanno.panTEannot.gff3.cleaned.gff")
    log:
        "logs/{sample}_postprocess_final.log"
    threads: config["resources"]["threads"]
    params:
        directory=W("single_genome_EDTA/{sample}")
    shell:
        """
        source activate tools
        scripts_for_snakepipe/postprocess_panTE_annotation.sh {input.bed} {input.gff} {input.trashgff} > {log} 2>&1 
        """