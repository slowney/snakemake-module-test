from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider()

rule all:
    input:
        OutputBam=expand("Aligned_BAMs/{sample}.bam", sample=config["Samples"]),


rule Align:
    input:
        # InputFQ1=rules.ExtractUMI.output.OutputFQ1,
        # InputFQ2=rules.ExtractUMI.output.OutputFQ2,
        InputFQ1=lambda wildcards: S3.remote(config["Samples"][wildcards.sample]["R1_Path"],keep_local=True),
        InputFQ2=lambda wildcards: S3.remote(config["Samples"][wildcards.sample]["R2_Path"]),
        # downloadComplete=rules.LocalizeGenome.output.downloadComplete,
        ref=S3.remote(multiext("s3://agc-027151828055-us-west-2/test_data/hg38_index/hg38.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"))
    output:
        OutputBam="Aligned_BAMs/{sample}.bam"
    params:
        Sample="{sample}",
        idx=lambda w, input: os.path.splitext(input.ref[0])[0]
    threads: 12
    resources:
        mem_mb=72000,
        disk_mb=100000
    conda: "envs/simpleAGC.yaml"
    shell:
        """
        echo INSTANCE: 
        curl http://169.254.169.254/latest/meta-data/instance-type
        
        bwa mem {params.idx} {input.InputFQ1} {input.InputFQ2} | samtools view -@ {threads} -Sb - > {output.OutputBam}
        """
