import pandas as pd

configfile: 'config.yaml'

localrules: all, link_fast5, fastq_to_fasta, merge_fasta

fast5_files = pd.read_table('doc/raw_reads.tsv').set_index('filename', drop=False)

def aggregate_fasta_input(wildcards):
    checkpoint_output = checkpoints.basecalling.get(**wildcards).output[0]
    print(checkpoint_output)
    gwc = glob_wildcards(os.path.join(checkpoint_output, 'fastq_runid_{run_id}_{batch_counter}_{part}.fastq'))
    print(gwc)
    return expand('fasta/{run}/fastq_runid_{run_id}_{batch_counter}_{part}.fa', zip,
                  run=[os.path.basename(checkpoint_output) for _ in gwc.run_id],
                  run_id=gwc.run_id,
                  batch_counter=gwc.batch_counter,
                  part=gwc.part)

rule all:
    input: expand('fasta/{run}.fa', run=set(fast5_files.run))

rule merge_run_fasta:
    input: aggregate_fasta_input
    output: 'fasta/{run}.fa'
    shell: 'cat {input} > {output}'

rule fastq_to_fasta:
    input: 'guppy_output/{run}/fastq_runid_{run_id}_{batch_counter}.fastq'
    output: 'fasta/{run}/fastq_runid_{run_id}_{batch_counter}.fa'
    conda: 'environment.yaml'
    shell: 'seqtk seq -A {input} > {output}' 

checkpoint basecalling:
    input: lambda wildcards: \
        expand('fast5/{run}/{filename}', \
               run=wildcards.run, \
               filename=[x for x, y in zip(fast5_files.filename, fast5_files.run) if y == wildcards.run])
    output: directory('guppy_output/{run}')
    wildcard_constraints:
        run='run\d+'
    params:
        flowcell='FLO-MIN106',
        kit='SQK-RAD004',
        parallel_callers=8
    threads: 4
    shell: 'guppy_basecaller '
           '--flowcell {params.flowcell} '
           '--kit {params.kit} '
           '--num_callers {params.parallel_callers} '
           '--cpu_threads_per_caller {threads} '
           '--save_path {output} '
           '--input_path fast5/{wildcards.run}'

rule link_fast5:
    input: lambda wildcards: \
        '{dirname}/{{filename}}'.format(dirname=fast5_files.dirname[wildcards.filename])
    output: 'fast5/{run}/{filename}'
    shell: 'ln -s {input} fast5/{wildcards.run}/'
