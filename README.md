# 2025_LifeOutdoors_Biodiversity_Workshop
## 🧬 Demographic and Diversity Analyses

This section outlines how to estimate historical demography using **PSMC**, identify **runs of homozygosity (ROH)**, and compute **genome-wide heterozygosity** from genomic data.

---

### 📊 PSMC (Pairwise Sequentially Markovian Coalescent)

#### 🛠️ Installation

```bash
# Create environment for running PSMC
conda create -n psmc -c bioconda psmc
conda activate psmc
conda install -c conda-forge gnuplot texlive-core ghostscript

# Create environment for plotting PSMC output
conda create -n psmc_plot
conda activate psmc_plot
conda install -c conda-forge/label/cf201901 gnuplot
conda install -c conda-forge texlive-core ghostscript
```

#### ▶️ Run PSMC

In your job script:

```bash
source ~/.bashrc
conda activate psmc
psmc -N25 -t15 -r5 -p "2+2+25*2+4+6" -o output.psmc output.psmcfa
#For bootstraping
splitfa output.psmcfa > split_output_psmcfa

seq 1 100 | xargs -P40 -I{} \
  psmc -N25 -t15 -r5 -b -p "2+2+25*2+4+6" \
       -o round-{}.psmc split_output_psmcfa

cat round-*.psmc > combined.psmc
```

#### 📈 Plot the PSMC results (in command line)

```bash
conda activate psmc_plot
perl psmc_plot.pl -R -p -u 2e-9 -g 1.0 psmc_output combined.psmc
```

> ⚠️ **Note**: Change `-g` to the generation time for your species. You can also use psmc_color_plot.pl (in each species directory) to change the plot colors.

---

### 🧬 Runs of Homozygosity (ROH)

#### 🛠️ Installation

```bash
# For ROH and genome assembly statistics
conda create -n bcftools -c bioconda bcftools assembly-stats

# For ROH using PLINK
conda create -n plink -c bioconda plink
```

#### ▶️ Run ROH analysis with bcftools (in a job script!)

```bash
bcftools roh \
  -G30 \
  --threads 24 \
  --AF-dflt 0.4 \
  filtered.vcf > sample.roh
```

#### ▶️ Run ROH analysis with PLINK (in terminatl

```bash
plink --vcf filtered.vcf --make-bed --double-id --allow-extra-chr --out sample

plink --bfile sample --homozyg --homozyg-snp 10 --homozyg-kb 20 --homozyg-density 1000 --homozyg-gap 10000 --homozyg-window-snp 10 --homozyg-window-het 2 --homozyg-window-missing 10 --homozyg-window-threshold 0.005 --allow-extra-chr --out sample_roh_relaxed
```

#### 📏 Estimate proportion of genome in ROH (FROH; in terminal )

```bash
# Total genome size
total_assembly=$(awk '/^>/ {if (seqlen) {sum += seqlen}; seqlen = 0; next} {seqlen += length} END {sum += seqlen; print sum}' ncbi_cleaned.fasta)

# Total ROH length
total_roh=$(awk '{sum += $9 * 1000} END {print sum}' sample_roh_relaxed.hom)

# Compute FROH
echo "scale=6; $total_roh / $total_assembly" | bc
```

This yields the inbreeding coefficient **FROH** as the fraction of the genome in homozygous tracts.

---

### 🧬 Genome-wide Heterozygosity

#### ▶️ Estimate with ANGSD and realSFS

First, run ANGSD to generate the `.saf.idx` file (not shown here), then run (in a job script):

```bash
realSFS angsdput.saf.idx > est.ml
```

To calculate heterozygosity (in terminal):

```bash
awk '{
  sum = $1 + $2 + $3;
  het = $2 / sum * 100;
  printf "heterozygosity = %.4f%%\n", het
}' est.ml
```

This gives the estimated genome-wide heterozygosity for one diploid individual.
