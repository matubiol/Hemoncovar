<div align="center">
  <h1> Hemoncovar </h1>
  <p><strong>A Comprehensive High-Throughput Snakemake Pipeline for Haemato-Oncology NGS Panels</strong></p>

  <p>
    <img src="https://img.shields.io/badge/Workflow-Snakemake-emerald?style=flat-square&logo=snakemake" alt="Snakemake">
    <img src="https://img.shields.io/badge/Language-Python%20%7C%20R%20%7C%20Bash-blue?style=flat-square&logo=python" alt="Languages">
    <img src="https://img.shields.io/badge/Genome-GRCh37%20%2F%20hg19-orange?style=flat-square" alt="Genome Assembly">
    <img src="https://img.shields.io/badge/License-MIT-yellow?style=flat-square" alt="License">
    <img src="https://img.shields.io/badge/Platform-HPC%20%7C%20Slurm-lightgrey?style=flat-square" alt="Platform">
  </p>
</div>

<hr />

<h2> Overview</h2>
<p>
  <code>Hemoncovar</code> is a high-throughput bioinformatic workflow for the analysis of targeted Next-Generation Sequencing (NGS) data in haemato-oncology (AML, CML, MPN, ALL, MDS). It was designed around the HealthInCode <em>Haematology OncoKit HS2</em> (89 genes, 26 fusion genes with intronic baits, genome-wide CNV backbone, dual UMIs; DNA from peripheral blood or bone marrow), but works with any capture panel given its BED file. Tumour-only somatic calling in blood and bone marrow brings specific challenges: no matched normal, subclonal mutations, difficult variants (FLT3-ITD, CALR, CEBPA, homopolymer indels) and recurrent fusions/rearrangements.
</p>
<p>
  This pipeline automates the entire process from reference infrastructure setup to fully annotated, clinical-grade variant reports by integrating state-of-the-art tools for small variants (SNVs/Indels), structural variations (SVs), and copy number alterations (CNVs).
</p>

<h2> Key Features & Modular Architecture</h2>

<table>
  <thead>
    <tr>
      <th>Module</th>
      <th>Tools Integrated</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><strong>0. Reference Preparation</strong></td>
      <td></td>
      <td>Harmonizes the chromosome naming convention of the intervals bed file with that of the reference fasta (e.g. "chr1" vs "1").</td>
    </tr>
    <tr>
      <td><strong>1. Quality Control </strong></td>
      <td><code>FastQC</code>, <code>Fastp</code>, <code>MultiQC</code></td>
      <td>Quality control on raw DNA-seq reads, including FastQC for quality assessment, Fastp for trimming and filtering, and MultiQC for aggregating QC reports.</td>
    </tr>
    <tr>
      <td><strong>2. Alignment </strong></td>
      <td><code>BWA-MEM2</code>, <code>GATK</code>, <code>MultiQC</code></td>
      <td>Aligns the processed reads to the reference genome using BWA-MEM2 and performs subsequent processing with GATK.</td>
    </tr>
    <tr>
      <td><strong>3. Small Somatic Variants</strong></td>
      <td><code>GATK Mutect2</code>, <code>VarDict</code>, <code>SnpEff</code></td>
      <td>Somatic variant calling using GATK Mutect2 and VarDict, followed by annotation with SnpEff and ClinVar. It includes configurable parameters for Mutect2 to allow users to adjust sensitivity and specificity based on their needs.</td>
    </tr>
    <tr>
      <td><strong>4. Structural Variants & Fusions</strong></td>
      <td><code>Manta</code>, <code>Lumpy</code></td>
      <td>Detection of genomic rearrangements, translocations, and pathologically relevant gene fusions, followed by gene annotation and clinical annotation using the CIViC database.</td>
    </tr>
    <tr>
      <td><strong>5. Copy Number Profiling</strong></td>
      <td><code>CNVkit</code>, <code>GATK</code></td>
      <td>Copy number profiling with CNVkit and GATK CNV against a reusable panel of normals (with leave-one-out QC), heterozygous-SNP B-allele frequencies, allelic imbalance / copy-neutral LOH and a per-sample report with arm-level events in cytogenetic notation.</td>
    </tr>
  </tbody>
</table>

<hr />

<h2> Prerequisites & Software</h2>
<p>All software is deployed automatically by Snakemake through conda; no environment modules or site-specific containers are needed.</p>
<ul>
  <li><strong>Snakemake wrappers</strong> (release set in <code>wrapper_version</code>, default <code>v9.18.0</code>) for fastp, FastQC, MultiQC, BWA-MEM2, samtools, Picard, fgbio (GroupReadsByUmi, CallMolecularConsensusReads, FilterConsensusReads), GATK (Mutect2 and filtering, BQSR, CNV tools), VarDict, bcftools, SnpEff/SnpSift and CNVkit.</li>
  <li><strong>Pinned conda environments</strong> in <code>envs/</code> for steps without a suitable wrapper, using the same tool versions as the wrapper release: fgbio FastqToBam/CorrectUmis/ZipperBams/CallDuplexConsensusReads and the streamed UMI re-alignments, Manta in tumour-only mode (the wrapper only supports <code>--normalBam</code>), Lumpy, GATK PreprocessIntervals/CreateReadCountPanelOfNormals/plots, CNVkit scatter, SnpSift extractFields and the Python/R scripts in <code>Pre_processing/</code> (BED harmonization, UMI length normalization) and <code>Post_processing/</code> (reports and annotation).</li>
</ul>

<details>
  <summary><b>Click to expand installation steps</b></summary>

  <h3>1. Clone the Repository</h3>
  <pre><code>git clone https://github.com/matubiol/Hemoncovar.git
cd Hemoncovar</code></pre>

  <h3>2. Install Snakemake</h3>
  <pre><code>mamba create -c conda-forge -c bioconda -n snakemake "snakemake>=8" snakemake-executor-plugin-slurm pandas
mamba activate snakemake</code></pre>

  <h3>3. (HPC without internet on compute nodes) Pre-build the environments</h3>
  <pre><code>snakemake -s Snakefile.smk --configfile config.yml --sdm conda --conda-create-envs-only --conda-prefix /shared/snakemake-conda</code></pre>
  <p>If the cluster has no access to GitHub, clone <code>https://github.com/snakemake/snakemake-wrappers</code>, check out the tag set in <code>wrapper_version</code> and add <code>--wrapper-prefix file:///path/to/snakemake-wrappers/</code>.</p>
</details>

<hr />

<h2> Reference & Data Preparation</h2>
<p>
  Before running the Snakemake workflow, you must download and index the GRCh37 human reference genome, resource bundles (dbSNP, Mills, gold standard indels), variant calling files (gnomAD panels), and clinical databases (CIViC, ClinVar).
</p>
<p>
  A comprehensive reference builder script is included in the repository. To set up the infrastructure:
</p>

<ol>
  <li>Open <code>src/BuildRef_GRCh37.sh</code> and configure your desired output directory (<code>OUTDIR</code>).</li>
  <li>Run the script using an environment with access to cluster modules:</li>
</ol>

<pre><code>bash src/BuildRef_GRCh37.sh</code></pre>

<p>
  This script will automatically index the reference genome, build the sequence dictionaries, intersect and format the GTF annotations with gene names instead of Ensembl IDs, filter gnomAD to your targeted regions via <code>bcftools</code>, pre-process GATK CNV target intervals, and cache the <strong>CIViC Database</strong>.
</p>

<hr />

<h2> Configuration</h2>
<p>Edit the <code>config.yml</code> file to point to the resource directory generated by your build script. Besides the files listed there, the following resources are required:</p>
<ul>
  <li><code>refs.common_snps</code>: common biallelic SNPs (e.g. <code>small_exac_common_3.vcf.gz</code>), ideally restricted to the panel targets. Used by <code>GetPileupSummaries</code> (Mutect2 contamination), <code>CollectAllelicCounts</code> (GATK CNV) and the heterozygous SNP genotyping for CNVkit (bcftools).</li>
  <li><code>refs.genome</code>: the reference fasta with its <code>.fai</code>, <code>.dict</code> and bwa-mem2 index (<code>.0123</code>, <code>.amb</code>, <code>.ann</code>, <code>.bwt.2bit.64</code>, <code>.pac</code>) next to it.</li>
  <li><code>refs.snpeff_data</code>: SnpEff data directory. The database <code>snpeff_db</code> is downloaded there automatically if it is missing.</li>
  <li><code>refs.gtf</code>: Ensembl GTF with <code>gene_name</code> and <code>gene_biotype</code> attributes (protein-coding genes are used for SV and CNV annotation).</li>
  <li><code>fusion_gene_panel</code> (optional): text file with one gene symbol per line (e.g. the 26 HS2 fusion genes) to restrict the SV/fusion report.</li>
</ul>

<h3>Sample sheet</h3>
<p>Samples are defined in a tab-separated sample sheet (<code>samplesheet</code> in <code>config.yml</code>, see <code>config/samples.tsv</code>):</p>
<table>
  <thead><tr><th>Column</th><th>Required</th><th>Description</th></tr></thead>
  <tbody>
    <tr><td><code>sample</code></td><td>yes</td><td>Sample identifier (underscores allowed).</td></tr>
    <tr><td><code>fastq_1</code>, <code>fastq_2</code></td><td>yes</td><td>R1/R2 FASTQs, absolute or relative to <code>raw_dir</code>. Several lanes are given comma-separated and are concatenated.</td></tr>
    <tr><td><code>type</code></td><td>no</td><td><code>tumor</code> (default) or <code>normal</code>. Normal samples (healthy donors processed with the same panel) build the CNV panels of normals and are not reported.</td></tr>
    <tr><td><code>sex</code></td><td>no</td><td><code>male</code> / <code>female</code>. Used for chrX/chrY in CNV calling.</td></tr>
    <tr><td><code>purity</code></td><td>no</td><td>Tumour fraction in (0, 1], e.g. blast percentage / 100. Used by CNVkit; if empty, no purity rescaling is applied.</td></tr>
  </tbody>
</table>
<p>If <code>samplesheet</code> is left empty, FASTQs in <code>raw_dir</code> are auto-discovered (Illumina naming <code>&lt;sample&gt;_S#_L00#_R1_001.fastq.gz</code> or <code>&lt;sample&gt;_R1.fastq.gz</code>), all lanes are kept, and the samples listed in <code>pon_samples</code> are treated as normals.</p>

<h3>CNV panel of normals</h3>
<p>CNV calling compares every tumour sample with a panel of normals (PoN) built from healthy donors sequenced with the same panel and protocol (at least 2, 10 or more recommended). Settings are in the <code>cnv:</code> block of <code>config.yml</code>:</p>
<ul>
  <li><code>pon_mode: build</code>: the PoN is built from the samples with <code>type=normal</code> and written as a reusable bundle to <code>results/03_VariantCalling/03_CNVs/00_PanelOfNormals</code> (CNVkit and GATK files plus <code>pon_manifest.tsv</code> with normals, sexes, versions and md5 checksums). To build only the PoN: <code>snakemake ... build_pon</code>.</li>
  <li><code>pon_mode: prebuilt</code>: reuse a bundle (<code>pon_dir</code>). The pipeline stops if the bundle was built with a different intervals BED. Rebuild it when the panel, kit lot or protocol changes.</li>
  <li><code>pon_qc</code> (build mode, 3 or more normals): each normal is analysed against a PoN built without it. <code>pon_qc.tsv</code> reports noise, false CNV calls, declared vs inferred sex and possible mosaic loss of chrY, with a keep/review recommendation. Remove flagged normals and rebuild.</li>
  <li>GATK excludes chrX/chrY by default (<code>gatk_exclude_sex_chroms</code>), as the PoN usually mixes sexes; sex chromosomes are reported by CNVkit, which is sex-aware when the sample sheet provides the sex.</li>
  <li>Check that the panel BED includes the genome-wide backbone probes: without them, arm-level CNV and LOH can only be assessed around the panel genes.</li>
</ul>

<h3>Main outputs</h3>
<ul>
  <li><code>03_VariantCalling/01_SmallVariants/&lt;caller&gt;/02_Annotated/&lt;sample&gt;_report.tsv</code>: coordinates, REF/ALT, gene, transcript, HGVS, effect, depth, alt depth, VAF, caller FILTER, ClinVar and a preliminary class (ClinVar or predicted effect; to be reviewed, not an ACMG/AMP classification). Only PASS calls by default (<code>report_keep_filtered</code>).</li>
  <li><code>03_VariantCalling/02_StructuralVariants/&lt;caller&gt;/02_Annotated/&lt;sample&gt;_ann_civic.xlsx</code>: SVs with genes at both breakpoints, fusion candidates (including intrachromosomal DEL/INV/DUP fusions such as STIL::TAL1 or CBFB::MYH11), genes spanned by deletions/duplications and CIViC diseases.</li>
  <li><code>03_VariantCalling/03_CNVs/Report/&lt;sample&gt;_CNV_report.xlsx</code>: CNV report combining CNVkit and GATK (sheets Summary, Arm_level, Segments, Genes, Allelic_imbalance, Parameters).</li>
  <li><code>03_VariantCalling/03_CNVs/&lt;caller&gt;/.../&lt;sample&gt;.ann.tsv</code>: all segments of each caller with log2 ratio, call, genes, cytoband and arm; CNVkit scatter/diagram and GATK plots.</li>
</ul>

<hr />

<h2> Running the Pipeline</h2>

<h3>Execution</h3>
<p>To run the pipeline on a Slurm cluster:</p>
<pre><code>snakemake -s Snakefile.smk --configfile config.yml --profile slurm --sdm conda --conda-prefix /shared/snakemake-conda</code></pre>
<p>Slurm profiles can be configured in <code>config/slurm/</code> to specify partition, account, memory, and other cluster-specific parameters. Other options can be specified as needed. Find more information in the <a href="https://snakemake.readthedocs.io/en/stable/executing/cli.html#defining-global-profiles">Snakemake documentation</a>.</p>

<hr />

<h2> License</h2>
<p>This project is licensed under the <strong>MIT License</strong> - see the <a href="LICENSE">LICENSE</a> file for details. It grants free permission to use, modify, distribute, and commercially exploit the software, provided the original copyright notice is preserved.</p>

<hr />

<div align="center">
  <p>Developed for Hemato-Oncology Research Groups. For issues, bug reports, or feature requests, please open a GitHub Issue.</p>
</div>
