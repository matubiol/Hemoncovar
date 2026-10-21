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
  <code>Hemoncovar</code> is a high-throughput bioinformatic workflow for the analysis of targeted Next-Generation Sequencing (NGS) data in haemato-oncology (AML, CML, MPN, ALL, MDS). It was designed around the HealthinCode <em>Haematology OncoKit HS2</em> (89 genes, 26 fusion genes with intronic baits, genome-wide CNV backbone, dual UMIs; DNA from peripheral blood or bone marrow), but works with any capture panel given its BED file. Tumour-only somatic calling in blood and bone marrow brings specific challenges: no matched normal, subclonal mutations, difficult variants (FLT3-ITD, CALR, CEBPA, homopolymer indels) and recurrent fusions/rearrangements.
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

  <h3>3. (Optional) Pre-build the conda environments</h3>
<p>
  Snakemake creates the environments automatically on the first run.
  To build them beforehand (e.g. to share them between users), run:
</p>
<pre><code>snakemake -s Snakefile.smk --configfile config.yml --sdm conda --conda-create-envs-only --conda-prefix /shared/snakemake-conda</code></pre>
<p>
  Use the same <code>--conda-prefix</code> in later runs.
</p>
<p>
  If compute nodes have no access to GitHub, clone
  <a href="https://github.com/snakemake/snakemake-wrappers">https://github.com/snakemake/snakemake-wrappers</a>,
  check out the tag set in <code>wrapper_version</code> and add
  <code>--wrapper-prefix file:///path/to/snakemake-wrappers/</code>.
</p>
</details>

<hr />

<h2> Reference & resources </h2>
<p>
  Before running the Snakemake workflow, you must download and index the human reference genome, resource bundles (dbSNP, Mills, gold standard indels), variant calling files (gnomAD panels), and clinical databases (CIViC, ClinVar). As an example, the following sources were used for a GRCh37/b37 setup; equivalent files for other assemblies (e.g. GRCh38) can be obtained from the same providers:
</p>
<ul>
  <li><strong>Human reference genome (GRCh37, Ensembl release 75):</strong> <a href="https://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz">Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz</a></li>
  <li><strong>Gene annotation (GTF, Ensembl GRCh37 release 87):</strong> <a href="https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz">Homo_sapiens.GRCh37.87.chr.gtf.gz</a></li>
  <li><strong>Base recalibration resource bundle</strong> (<a href="https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references;tab=objects?prefix=&forceOnObjectsSortingFiltering=false">Broad Institute GATK resource bundle, b37</a>):
    <ul>
      <li>dbSNP: <a href="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/dbsnp_138.b37.vcf.gz">dbsnp_138.b37.vcf.gz</a> (<a href="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/dbsnp_135.b37.vcf.gz.tbi">.tbi</a>)</li>
      <li>Mills and 1000G gold standard indels: <a href="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Mills_and_1000G_gold_standard.indels.b37.vcf.gz">Mills_and_1000G_gold_standard.indels.b37.vcf.gz</a> (<a href="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Mills_and_1000G_gold_standard.indels.b37.vcf.gz.tbi">.tbi</a>)</li>
      <li>1000G phase 1 high-confidence SNPs: <a href="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/1000G_phase1.snps.high_confidence.b37.vcf.gz">1000G_phase1.snps.high_confidence.b37.vcf.gz</a> (<a href="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/1000G_phase1.snps.high_confidence.b37.vcf.gz.tbi">.tbi</a>)</li>
    </ul>
  </li>
  <li><strong>Variant calling resources</strong> (<a href="https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38;tab=objects?prefix=&forceOnObjectsSortingFiltering=false">GATK somatic best practices, b37</a>):
    <ul>
      <li>gnomAD allele-frequency-only sites: <a href="https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf">af-only-gnomad.raw.sites.vcf</a> (<a href="https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf.idx">.idx</a>)</li>
      <li>Mutect2 exome panel of normals: <a href="https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf">Mutect2-exome-panel.vcf</a> (<a href="https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf.idx">.idx</a>)</li>
      <li>Common germline variants (ExAC): <a href="https://storage.googleapis.com/gatk-best-practices/somatic-b37/small_exac_common_3.vcf">small_exac_common_3.vcf</a> (<a href="https://storage.googleapis.com/gatk-best-practices/somatic-b37/small_exac_common_3.vcf.idx">.idx</a>)</li>
    </ul>
  </li>
  <li><strong>CIViC</strong> clinical evidence summaries (01-Apr-2026 release): <a href="https://civicdb.org/downloads/01-Apr-2026/01-Apr-2026-AcceptedAndSubmittedClinicalEvidenceSummaries.tsv">01-Apr-2026-AcceptedAndSubmittedClinicalEvidenceSummaries.tsv</a></li>
  <li><strong>ClinVar</strong> (GRCh37): <a href="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz">clinvar.vcf.gz</a> (<a href="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi">.tbi</a>)</li>
  <li><strong>Cytogenetic bands</strong> (UCSC hg19): <a href="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz">cytoBand.txt.gz</a></li>
</ul>

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
<pre><code>snakemake -s Snakefile.smk --configfile config.yml --profile slurm</code></pre>
<p>An example Snakemake profile is provided in <code>config/config.yaml</code>. It sets the Slurm executor, the number of parallel jobs, the software deployment methods (Conda, environment modules and Apptainer) and the default cluster resources. Before running the pipeline, adapt the following fields to your cluster:</p>
<ul>
  <li><code>slurm_account</code> and <code>slurm_partition</code>: your Slurm account and partition.</li>
  <li><code>mem_mb</code>: default memory per job (in MB).</li>
  <li><code>jobs</code>: maximum number of jobs submitted simultaneously.</li>
  <li><code>conda-prefix</code>: shared directory where Conda environments will be created and cached.</li>
  <li><code>apptainer-args</code>: paths to bind inside the containers.</li>
</ul>
<p>Any option defined in the profile can be overridden from the command line (e.g. <code>--jobs 20</code>). Other options can be added as needed. Find more information in the <a href="https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles">Snakemake documentation</a>.</p>

<hr />

<h2> License</h2>
<p>This project is licensed under the <strong>MIT License</strong> - see the <a href="LICENSE">LICENSE</a> file for details. It grants free permission to use, modify, distribute, and commercially exploit the software, provided the original copyright notice is preserved.</p>

<hr />

<div align="center">
  <p>Developed for Hemato-Oncology Research Groups. For issues, bug reports, or feature requests, please open a GitHub Issue.</p>
</div>
