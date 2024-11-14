import os
import glob

# constants
MAFIA_DIR = "/data/gpfs/projects/punim0614/occheng/m6A-QTL/mAFiA"

BASE_DIR = "/data/gpfs/projects/punim0614/shared/shared_data/external_public/PRJEB76585_directRNA/media/nsl114/Nanopore/ONT_RNApolyA_24112020/mRNAONT24112020"
FASTQ_DIR = BASE_DIR + "/{sample}/{seq}/fastq_pass"

# use transcriptome fa
REF = "/data/gpfs/projects/punim0614/shared/shared_data/external_public/refgenome_annotation/GENCODE/GRCh37.p13/gencode.v19.pc_transcripts.fa"

def get_sample_dir(wc, suffix = ""):
  ff = glob.glob(BASE_DIR + "/%s/*%s" % (wc.sample, suffix))
  if len(ff) != 1:
    raise Exception('Expected exactly 1 file for kind "%s" - got %s' % (wc.sample, ff))
  
  return ff[0]

def get_fastq_dir(wc):
  return get_sample_dir(wc, suffix = "/fastq_pass")
  

CHROMOSOMES = list(range(1, 23)) + ["X", "Y"]

# This rule runs everything needed for all the samples in BASE_DIR
(SAMPLES,SEQS) = glob_wildcards(FASTQ_DIR)
#SAMPLES,SEQS = (SAMPLES[0], SEQS[0])
rule all:
  input:
    expand("{sample}/m6anet/data.indiv_proba.csv", sample=SAMPLES),
    #expand("{sample}/mAFiA/mapped/mapped.bam", sample=SAMPLES),
    #expand("{sample}/mAFiA/predicted/chr{chr}", sample=SAMPLES, chr=CHROMOSOMES)
    # expand("{sample}/cheui/site_level_m6A_predictions.txt", sample=SAMPLES)


rule map_and_sort:
  input: get_fastq_dir
  output:
    bam = "{sample}/mapped/pass_alignment.bam",
    sorted_bam = "{sample}/mapped/pass_alignment.sorted.bam",
    index = "{sample}/mapped/pass_alignment.sorted.bam.bai"
  envmodules:
    "GCCcore/11.3.0",
    "GCC/11.3.0",
    "minimap2/2.26",
    "SAMtools/1.16.1",
  threads: 16
  resources:
    mem_mb=32000,
    time="1-00:00:00",
    runtime="1d"
  shell:
    """
    minimap2 -t {threads} -ax splice -k14 -t -uf {REF} {input}/*.fastq \
        | samtools view -Sb - > {output.bam}
    samtools sort -@ {threads} -O bam -o {output.sorted_bam} {output.bam}
    samtools index {output.sorted_bam}
    """
    
rule merge_basecalled_into_single_fastq:
  input: get_fastq_dir
  output:
    "{sample}/nanopolish/merged_reads.fastq"
  shell:
    "cat {input}/*.fastq > {output}"
    
rule nanopolish_index:
  input:
    fast5 = lambda wc: get_sample_dir(wc, "/fast5_pass/"),
    fastq = "{sample}/nanopolish/merged_reads.fastq"
  output:
    "{sample}/nanopolish/merged_reads.fastq.index",
    "{sample}/nanopolish/merged_reads.fastq.index.fai",
    "{sample}/nanopolish/merged_reads.fastq.index.gzi"
    #"{sample}/nanopolish/merged_reads.fastq.index.readdb"
  envmodules:
    "GCC/11.3.0",
    "OpenMPI/4.1.4",
    "nanopolish/0.13.3"
  resources:
    mem_mb=16000,
    time="6:00:00",
    runtime="6h"
  shell:
    "nanopolish index -d {input.fast5} {input.fastq}"
    
summary_fn = lambda wc: glob.glob(get_sample_dir(wc, "/sequencing_summary_*.txt"))
    
rule nanopolish_eventalign_m6anet:
  input:
    reads = "{sample}/nanopolish/merged_reads.fastq",
    bam = "{sample}/mapped/pass_alignment.sorted.bam",
    bai = "{sample}/mapped/pass_alignment.sorted.bam.bai",
    index = "{sample}/nanopolish/merged_reads.fastq.index",
    summary = summary_fn
  output:
    "{sample}/nanopolish/pass_alignment.eventalign.txt"
  threads: 96
  resources:
    mem_mb=32000,
    time="1-12:00:00",
    runtime="36h"
  envmodules:
    "GCC/11.3.0",
    "OpenMPI/4.1.4",
    "nanopolish/0.13.3"
  shell:
    "nanopolish eventalign "
    "--reads {input.reads} --bam {input.bam} "
    "--genome {REF} --scale-events "
    "--threads={threads} --signal-index "
    "--summary {input.summary} > {output}"
    
# TODO: are these the same??
rule nanopolish_eventalign_cheui:
  input:
    reads = "{sample}/nanopolish/merged_reads.fastq",
    bam = "{sample}/mapped/pass_alignment.sorted.bam",
    bai = "{sample}/mapped/pass_alignment.sorted.bam.bai",
    index = "{sample}/nanopolish/merged_reads.fastq.index",
    summary = summary_fn
  output:
    "{sample}/nanopolish/pass_alignment.eventalign.with_names.txt"
  threads: 64
  resources:
    mem_mb=32000,
    time="1-12:00:00",
    runtime="36h"
  envmodules:
    "GCC/11.3.0",
    "OpenMPI/4.1.4",
    "nanopolish/0.13.3"
  shell:
    "nanopolish eventalign -t {threads} "
    "--reads {input.reads} "
    "--bam {input.bam} "
    "--genome {REF} "
    "--scale-events --signal-index --samples --print-read-names > {output}"

#
#
# M6ANET
#
#

rule m6anet_dataprep:
  input:
    "{sample}/nanopolish/pass_alignment.eventalign.txt"
  output:
    "{sample}/m6anet/data.json",
    "{sample}/m6anet/data.log",
    "{sample}/m6anet/data.info",
    "{sample}/m6anet/eventalign.index"
  threads: 48
  resources:
    mem_mb=192000,
    time="24:00:00",
    runtime="24h"
  conda: "env/m6anet.yaml"
  shell:
    "echo Dataprep && m6anet dataprep --eventalign {input} "
    "--out_dir {wildcards.sample}/m6anet/ --n_processes {threads}"

rule m6anet_inference:
  input:
    "{sample}/m6anet/data.json"
  output:
    "{sample}/m6anet/data.indiv_proba.csv"
  conda: "env/m6anet.yaml"
  threads: 16
  resources:
    mem_mb=32000,
    time="4:00:00",
    runtime="4h"
  shell:
    "echo Starting && m6anet inference --input_dir {wildcards.sample}/m6anet "
    "--out_dir {wildcards.sample}/m6anet "
    "--n_processes {threads} --num_iterations 1000"
    
    
#
#
# CHEUI
#
#

rule cheui_preprocess:
  input:
    "{sample}/nanopolish/pass_alignment.eventalign.with_names.txt"
  output:
    dir = directory("{sample}/cheui/out_A_signals+IDs.p/"),
    file = "{sample}/cheui/out_A_signals+IDs.p/nanopolish_output_test_signals+IDS.p"
  threads: 64
  envmodules:
    "GCC/11.3.0"
  resources:
    mem_mb=64000,
    time="2-00:00:00",
    runtime="2d"
  shell:
    """
    # save input and output before cd
    INPUT=$(realpath "{input}")
    OUTPUT=$(realpath "{output.dir}")
    
    cd {CHEUI_DIR}/scripts/preprocessing_CPP
    ./CHEUI -i $INPUT -m ../../kmer_models/model_kmer.csv -n {threads} --m6A -o $OUTPUT
    """
    
rule cheui_readmodel:
  input: 
    "{sample}/cheui/out_A_signals+IDs.p/nanopolish_output_test_signals+IDS.p"
  output:
    "{sample}/cheui/read_level_m6A_predictions.txt"
  conda: "env/cheui.yaml"
  resources:
    partition="gpu-a100",
    gpus=1,
    mem_mb=32000,
    time="1-00:00:00",
    runtime="1d"
  shell:
    "python {CHEUI_DIR}/scripts/CHEUI_predict_model1.py "
    "-i {input} "
    "-m {CHEUI_DIR}/CHEUI_trained_models/CHEUI_m6A_model1.h5 "
    "-o {output} -l WT_rep1"
  
rule sort_read_level:
  input:
    "{sample}/cheui/read_level_m6A_predictions.txt"
  output:
    "{sample}/cheui/read_level_m6A_predictions.sorted.txt"
  threads: 4
  shell:
    "sort -k1 --parallel={threads} {input} > {output}"
    
rule cheui_sitemodel:
  input:
    "{sample}/cheui/read_level_m6A_predictions.sorted.txt"
  output:
    "{sample}/cheui/site_level_m6A_predictions.txt"
  conda: "env/cheui.yaml"
  resources:
    partition="gpu-a100",
    gpus=1,
    mem_mb=32000,
    time="1-00:00:00",
    runtime="1d"
  shell:
    "python3 {CHEUI_DIR}/scripts/CHEUI_predict_model2.py "
    "-i {input} "
    "-m {CHEUI_DIR}/CHEUI_trained_models/CHEUI_m6A_model2.h5 "
    "-o {output}"
    
#
#
# mAFiA
#
#

rule mafia_basecall:
  input:
    lambda wc: get_sample_dir(wc) + "/fast5_pass/"
  output:
    directory("{sample}/mAFiA/basecalled/")
  resources:
    slurm_extra="-p gpu-a100 --gres=gpu:1",
    mem_mb=32000,
    time="1-00:00:00",
    runtime="1d"
  shell:
    """
    source {MAFIA_DIR}/mafia-venv/bin/activate
    
    python3 {MAFIA_DIR}/RODAN/basecall.py \
      --fast5dir {input} \
      --model {MAFIA_DIR}/models/backbone.torch \
      --batchsize 4096 \
      --outdir {output}
    """
    
rule mafia_align:
  input:
    "{sample}/mAFiA/basecalled/rodan.fasta"
  output:
    bam = "{sample}/mAFiA/mapped/mapped.bam",
    index = "{sample}/mAFiA/mapped/mapped.bam.bai"
  envmodules:
    "GCCcore/11.3.0",
    "GCC/11.3.0",
    "minimap2/2.26",
    "SAMtools/1.16.1",
  threads: 16
  resources:
    mem_mb=32000,
    time="04:00:00",
    runtime="4h"
  shell:
    """
    minimap2 --secondary=no -ax splice -uf -k14 -t {threads} --cs {REF} {input} \
    | samtools view -bST {REF} -q50 - \
    | samtools sort - > {output.bam}
    
    samtools index {output.bam}
    """
    
BACKBONE = MAFIA_DIR + "/models/backbone.torch"
CLASSIFIERS = MAFIA_DIR + "/models/classifiers"
      
rule mafia_predict:
  input:
    bam = "{sample}/mAFiA/mapped/mapped.bam",
    fast5dir = lambda wc: get_sample_dir(wc) + "/fast5_pass/",
    chr = MAFIA_DIR + "/misc/m6A.GRCh37.chrchr{chr}.bed"
  output:
    directory("{sample}/mAFiA/predicted/chr{chr}")
  resources:
    slurm_extra="-p gpu-a100 --gres=gpu:1",
    mem_mb=128000,
    time="2-00:00:00",
    runtime="48h"
  shell:
    """
    source {MAFIA_DIR}/mafia-venv/bin/activate
    
    test_mAFiA \
    --bam_file {input.bam} \
    --fast5_dir {input.fast5dir} \
    --ref_file {REF} \
    --mod_file {input.chr} \
    --min_coverage 50 \
    --max_num_reads 1000 \
    --backbone_model_path {BACKBONE} \
    --classifier_model_dir {CLASSIFIERS} \
    --mod_prob_thresh 0.5 \
    --out_dir {output}
    """