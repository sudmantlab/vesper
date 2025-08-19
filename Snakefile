import glob
import os

def get_specimens():
    vcf_files = glob.glob("analysis/files/hg38_scaffolded/*.qc_all.vcf.gz")
    return [os.path.basename(f).replace(".qc_all.vcf.gz", "") for f in vcf_files]

rule all:
    input:
        expand("output/{specimen}/{specimen}.annotated.refined.vcf.gz", specimen=get_specimens())

rule pre_filter_vcf:
    input:
        "analysis/files/hg38_scaffolded/{specimen}.qc_all.vcf.gz"
    output:
        temp("output/{specimen}/{specimen}.pre.vcf.gz")
    params:
        svlen = 10000,
        max_cov = 60,
        min_cov = 20,
        max_freq = 0.1
    threads: 1
    shell:
        """
        regions=$(for i in {{1..22}} X Y; do echo -n "chr${{i}}_RagTag_hap1,chr${{i}}_RagTag_hap2,"; done | sed 's/,$//')
        bcftools view -r $regions {input} | \
        bcftools annotate --set-id '{wildcards.specimen}\_%ID' - | \
        bcftools filter -i 'FILTER="PASS"|FILTER="ALN_NM"' - | \
        bcftools filter -i 'INFO/SVTYPE="INS" & INFO/SVLEN <= {params.svlen} & {params.max_cov} > (DR + DV) & (DR + DV) > {params.min_cov} & (DV <= {params.max_freq} * (DR + DV))' - \
        -o {output} --write-index   
        """

rule annotate_vcf:
    input:
        "output/{specimen}/{specimen}.pre.vcf.gz"
    output:
        temp("output/{specimen}/{specimen}.pre.annotated.vcf.gz"),
        temp("output/{specimen}/{specimen}.pre.annotated.vcf.gz.tbi")
    threads: 8
    shell:
        """
        vesper annotate --vcf {input} --output-dir output/{wildcards.specimen} --threads {threads}
        """

rule filter_annotated_vcf:
    input: 
        "output/{specimen}/{specimen}.pre.annotated.vcf.gz"
    output:
        "output/{specimen}/{specimen}.annotated.vcf.gz"
    threads: 1
    shell:
        """
        bcftools filter -i "REPEATMASKER_RESULTS!='.'" {input} -o {output} --write-index
        """


rule refine_vcf:
    input:
        vcf = "output/{specimen}/{specimen}.annotated.vcf.gz",
        bam = "/Volumes/mnemosyne/spermSV/output/alignment/hg38_scaffolded/minimap2/standard/mapped/{specimen}.sorted.merged.bam"
    output:
        "output/{specimen}/{specimen}.annotated.refined.vcf.gz"
    threads: 4
    shell:
        """
        vesper refine --vcf {input.vcf} --bam {input.bam} --output-dir output/{wildcards.specimen} --threads {threads}
        """