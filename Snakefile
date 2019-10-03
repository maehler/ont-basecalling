import pandas as pd

configfile: 'config.yaml'

localrules: all, link_fast5

fast5_files = pd.read_table('doc/raw_reads.tsv')
linked_fast5_files = [os.path.join('fast5', os.path.basename(x)) for x in fast5_files.filename]

def aggregate_fasta_input(wildcards):
    checkpoint_output = checkpoints.basecalling.get(**wildcards).output[0]
    gwc = glob_wildcards(os.path.join(checkpoint_output, 'fastq_runid_{run_id}_{batch_counter}.fastq'))
    return expand('fasta/all/fastq_runid_{run_id}_{batch_counter}.fa', **gwc)

rule all:
    input: 'fasta/aspen.fa'

rule merge_fasta:
    input: aggregate_fasta_input
    output: 'fasta/aspen.fa'
    shell: 'cat {input} > {output}'

rule fastq_to_fasta:
    input: 'guppy_output/fastq_runid_{run_id}_{batch_counter}.fastq'
    output: 'fasta/all/fastq_runid_{run_id}_{batch_counter}.fa'
    conda: 'environment.yaml'
    shell: 'seqtk seq -A {input} > {output}' 

checkpoint basecalling:
    input: linked_fast5_files
    output: directory('guppy_output')
    params:
        flowcell='FLO-MIN106',
        kit='SQK-RAD004',
        parallel_callers=8,
        indir='fast5'
    threads: 4
    shell: 'guppy_basecaller '
           '--flowcell {params.flowcell} '
           '--kit {params.kit} '
           '--num_callers {params.parallel_callers} '
           '--cpu_threads_per_caller {threads} '
           '--save_path {output} '
           '--input_path {params.indir}'

rule link_fast5:
    input: fast5_files.filename
    output: linked_fast5_files
    shell: 'ln -s {input} fast5/'
