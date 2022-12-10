import os

configfile: 'config.yaml'

SAMPLEID = config['sample_ids']


bam_input = os.path.join(config['input_dir'], "{sampleid}.bam")
bam_ch22_input = os.path.join(config['input_dir'], "{sampleid}.ch22.bam")
   

rule all:
    input:
        expand(str(config['input_dir'] + "/{sampleid}.bam.bai"), sampleid = SAMPLEID),
        expand(str(config['output_dir'] + "/{sampleid}.wgs_metrics.txt"), sampleid = SAMPLEID),
        expand(str(config['output_dir'] + "/{sampleid}.gc_bias_metrics.txt"), sampleid = SAMPLEID),
        expand(str(config['output_dir'] + "/{sampleid}.gc_bias_metrics_chart.pdf"), sampleid = SAMPLEID),
        expand(str(config['output_dir'] + "/{sampleid}.gc_bias_summary.txt"), sampleid = SAMPLEID),
        expand(str(config['output_dir'] + "/{sampleid}.duplicate_metrics.txt"), sampleid = SAMPLEID),
        expand(str(config['output_dir'] + "/{sampleid}.insert_metrics.txt"), sampleid = SAMPLEID)
      


#index bam files for processing 
rule index_bam_files:
   input:
      bam_input
   output:
      index_output = os.path.join(config['input_dir'], "{sampleid}.bam.bai")
   shell:
    """ 
    samtools index {input}
    """

#pull out chromosome 22 from dataset 
rule pull_ch22:
    input: 
        bam_input
    output:
        os.path.join(config['input_dir'], "{sampleid}.ch22.sam")
    shell:
        """
        samtools view -h {input} 22 > {output}
        """

#compress ch22 sam file back down to bam file 
rule compress_sam:
    input:
        os.path.join(config['input_dir'], "{sampleid}.ch22.sam")
    output:
        os.path.join(config['input_dir'], "{sampleid}.ch22.bam")
    shell:
        """
        samtools view -bS {input} > {output}
        """

# Collect WGS Metrics 
rule wgs_metrics_single_end:
    input: 
        reference = config['reference_genome'],
        bam = bam_ch22_input
    output:
        os.path.join(config['output_dir'], "{sampleid}.wgs_metrics.txt")
    shell:
        """
        gatk CollectWgsMetrics \
            I={input.bam} \
            O={output} \
            R={input.reference}
        """

# Collect GC Bias Metrics 
rule collect_gc_bias_metrics:
    input: 
        bam = bam_ch22_input, 
        reference = config['reference_genome']
    output:
        metrics = os.path.join(config['output_dir'], "{sampleid}.gc_bias_metrics.txt"),
        chart = os.path.join(config['output_dir'], "{sampleid}.gc_bias_metrics_chart.pdf"),
        summary = os.path.join(config['output_dir'], "{sampleid}.gc_bias_summary.txt")
    shell:
        """
        gatk CollectGcBiasMetrics \
            I={input.bam} \
            O={output.metrics} \
            CHART={output.chart} \
            S= {output.summary} \
            R={input.reference}
        """

# Collect Duplicate Metrics
rule collect_duplicate_metrics:
    input:
        bam_ch22_input
    output:
        os.path.join(config['output_dir'], "{sampleid}.duplicate_metrics.txt")
    shell:
        """
        gatk CollectDuplicateMetrics \
            I= {input} \
            M= {output}
        """

# Collect Insert Size Metrics 
rule collect_insert_size_metrics:
    input:
        bam_ch22_input
    output:
        metrics = os.path.join(config['output_dir'], "{sampleid}.insert_metrics.txt"),
        histogram = os.path.join(config['output_dir'], "{sampleid}.insert_metrics_histo.pdf")
    shell:
        """
        gatk CollectInsertSizeMetrics \
          I= {input} \
          O= {output.metrics} \
          H= {output.histogram}  
        """