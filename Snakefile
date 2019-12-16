import pandas as pd

configfile: 'config.yaml'

shell.executable('bash')
shell.prefix('source {module_init}; set -euo pipefail; ' \
             .format(module_init=config['cluster']['module_init']))

localrules: all, link_fast5, fastq_to_fasta, merge_run_fasta

fast5_dfs = []
for fname in config['raw']:
    fast5_dfs.append(pd.read_table(fname).set_index('filename', drop=False))
fast5_files = pd.concat(fast5_dfs)

def aggregate_fasta_input(wildcards):
    checkpoint_output = checkpoints.gpu_basecalling.get(**wildcards).output[0]
    gwc = glob_wildcards(os.path.join(checkpoint_output, 'fastq_runid_{run_id}_{batch_counter}_{part}.fastq'))
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
    wildcard_constraints:
        run=r'[^\/]+'
    shell: 'cat {input} > {output}'

rule fastq_to_fasta:
    input: 'guppy_output/{run}/fastq_runid_{run_id}_{batch_counter}_{part}.fastq'
    output: 'fasta/{run}/fastq_runid_{run_id}_{batch_counter}_{part}.fa'
    wildcard_constraints:
        batch_counter=r'\d+',
        part=r'\d+'
    conda: 'environment.yaml'
    shell: 'seqtk seq -A {input} > {output}' 

checkpoint gpu_basecalling:
    input: lambda wildcards: \
        expand('fast5/{run}/{filename}', \
               run=wildcards.run, \
               filename=fast5_files.filename[fast5_files.run == wildcards.run])
    output: directory('guppy_output/gpu_{run}')
    wildcard_constraints:
        run=r'[^\/]+'
    params:
        flowcell=lambda wildcards: fast5_files.flowcell[fast5_files.run == wildcards.run][0],
        kit=lambda wildcards: fast5_files.chemistry[fast5_files.run == wildcards.run][0]
    shell:
        '''
        module load bioinfo-tools guppy
        guppy_basecaller \
            --flowcell {params.flowcell} \
            --kit {params.kit} \
            --gpu_runners_per_device 20 \
            --device "cuda:0" \
            --save_path {output} \
            --input_path fast5/{wildcards.run}
        '''

rule link_fast5:
    input: lambda wildcards: \
        '{dirname}/{{filename}}'.format(dirname=fast5_files.dirname[wildcards.filename])
    output: 'fast5/{run}/{filename}'
    shell: 'ln -s {input} fast5/{wildcards.run}/'
