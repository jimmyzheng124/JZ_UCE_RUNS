#Ultraconserved Element (UCE) Analyses for Malawi Cichlids
This github serves the purpose of running analyses on DNA ultra-conserved elements data for Malawi cichlids. You will find configuration files, pipeline scripts, and reconstructed trees.
Credit goes to Brant Faircloth for constructing the <a href="https://github.com/faircloth-lab/phyluce">phyluce pipeline</a>.


***Jimmy Zheng***  
*December 11, 2015*  

See [Final Matrices and Phylogenetic Analysis] (#phylogeny) to skip to the phylogenetic analysis pipeline, post-assembly, post-alignment, and post-trimming.

For the data preparation section immediately ahead, click on the following to navigate to specific steps: [processing raw reads] (#rawreads), [assembling data] (#assembly), [incorporating external data] (#outgroup), and [aligning sequences] (#alignment).

##<a name="preparation">Preparing UCE DATA for Phylogenetic Analysis</a>
We have generated UCE reads data for 90 species samples of cichlids, 60 from Lake Tanganyika and 30 from Lake Malawi. On average, we processed about 2 million reads for each species. You want between 1 to 2 million reads for good base coverage.

Here, we will analyze both datasets (see [Raw Reads] (#rawreads)) and eventually select only the Malawi cichlids for the final data matrices. The species are as follows:  
`*`*aristochromis christiae*  
`*`*astatotilapia burtoni [outgroup]*  
`*`*aulonocara stuartgranti*  
`*`*capidochromis eucinostomus*  
`*`*cheilotilapia euchilus*  
`*`*cheilotilapia rhodessi*  
`*`*copadichromis trimaculatus*  
`*`*cyathichromis obliquidens*  
`*`*docimodus evelynae*  
`*`*fossorichromis rostratus*  
`*`*genyochromis mento*  
`*`*hemitilapia oxyrhynchus*  
`*`*labeotropheus fullbornei*  
`*`*labeotropheus trewavase*  
`*`*labidochromis gigas*  
`*`*melanochromis auratus*  
`*`*melanochromis kaskazini*  
`*`*metriaclima greshakae*  
`*`*metriaclima patricki*  
`*`*nimbochromis polystigma*  
`*`*otopharynx heterodon*  
`*`*otopharynx lithobates*  
`*`*otopharynx picta*  
`*`*petrotilapia nigra*  
`*`*placidochromis electra*  
`*`*placidochromis milomo*  
`*`*psedotropheus crabro*  
`*`*pseudotropheus flavus*  
`*`*rhamphochromis longiceps*  
`*`*simochromis babaulti [outgroup]*  
`*`*stigmatochromis woodei*  
`*`*taeniolethrinops preorbitalis*  
`*`*tropheops microstoma*  
`*`*tyrannochromis nigriventer*

###<a name="rawreads">STEP 1: CREATE DIRECTORIES AND PROCESS RAW READS</a>
When doing UCE data analyses, you will want to keep all files as organized as possible. I first made a directory for all UCE analyses.  

```
mkdir JZ_UCE_RUNS; cd JZ_UCE_RUNS
mkdir raw-fastq; cd raw-fastq
unzip malawicichlids.zip
rm malawicichlids.zip
```
Next, I created a configuration file (.conf) for the illumiprocessor program, which will copy the raw reads into a different directory and remove contaminating adapter sequences and low-quality bases. This is where I match up each species to its respective adapter identifiers. See below. The following [file] (illumiprocessor.conf) includes both Tanganyika and Malawi cichlids.

```
[adapters]
#indexed adapters composed of a Y-yoked stub sequence and a unique identifier for both i7 and i5 ends
i7:GATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTG
i5:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT*GTGTAGATCTCGGTGGTCGCCGTATCATT

[tag sequences]
#the unique identifiers
i7-06_01:CACCTGTT
i7-07_01:ACCTTCTC
i7-07_02:ACACCAGT
i7-07_03:GAAGGAAG
i7-07_04:TCGAGTGA
...
i5-02_H:GCTTCGAA
i5-03_A:GTGGTGTT
i5-03_B:ACAGCTCA
i5-03_C:TTCCTGTG
i5-03_D:GGTTGTCA
...

[tag map]
#each species is tagged with a unique PAIR of adapters
aristochromis-christiae-ID25_S170:i7-08_01,i5-03_B
astatotilapia-burtoni-ID32_S177:i7-08_08,i5-03_B
aulonocara-stuartgranti-ID26_S171:i7-08_02,i5-03_B
bathybates-minor-ID90_S228:i7-08_07,i5-03_G
boulengerochromis-microlepis-ID79_S220:i7-08_08,i5-03_F
capidochromis-eucinostomus-ID12_S164:i7-07_12,i5-02_H
chalinochromis-brichardi-ID57_S202:i7-08_10,i5-03_D
cheilotilapia-euchilus-ID23_S152:i7-08_04,i5-03_H
...

[names]
#cleaning up sample names for the remainder of the analyses
aristochromis-christiae-ID25_S170:aristochromis-christiae-ID25
astatotilapia-burtoni-ID32_S177:astatotilapia-burtoni-ID32
aulonocara-stuartgranti-ID26_S171:aulonocara-stuartgranti-ID26
bathybates-minor-ID90_S228:bathybates-minor-ID90
boulengerochromis-microlepis-ID79_S220:boulengerochromis-microlepis-ID79
capidochromis-eucinostomus-ID12_S164:capidochromis-eucinostomus-ID12
chalinochromis-brichardi-ID57_S202:chalinochromis-brichardi-ID57
cheilotilapia-euchilus-ID23_S152:cheilotilapia-euchilus-ID23
... (continued)
```

Then I ran <a href="http://illumiprocessor.readthedocs.org/en/latest/usage.html">illumiprocessor</a> using the following code:

```
illumiprocessor \
    --input raw-fastq/ \
    --output dmux-merged-clean \
    --config illumiprocessor.conf \
    --cores 12
```


Your directory structure should now look like this:  

``` 
# Because directory structures from here on out become much more complicated, 
# I will not be displaying them frequently.
JZ_UCE_RUNS
├── dmux-merged-clean
│	 ├── aristochromis-christiae-ID25  
│	 │   ├── adapters.fasta
│	 │   ├── raw-reads
│	 │   ├── split-adapter-quality-trimmed
│	 │   │   ├── aristochromis-christiae-ID25-READ1.fastq.gz  
│	 │   │   ├── aristochromis-christiae-ID25-READ2.fastq.gz  
│	 │   │   └── aristochromis-christiae-ID25-READ-singleton.fastq.gz  
│	 │   └── stats  
│	 ├── harpagochromis-two-stripe-white-lip-ID88  
│	 ├── neolamprologus-bifasciatus-ID80  
│	 ├── otopharynx-lithobates-ID16  
│	 ├── astatotilapia-burtoni-ID32  
│	 └── ... (continued)
├── raw-fastq
│	 ├── aristochromis-christiae-ID25_S170_L999_R1_001.fastq.gz  
│	 ├── metriaclima-patricki-ID22_S151_L999_R1_001.fastq.gz  
│	 ├── aristochromis-christiae-ID25_S170_L999_R2_001.fastq.gz  
│	 ├── metriaclima-patricki-ID22_S151_L999_R2_001.fastq.gz  
│	 ├── astatotilapia-burtoni-ID32_S177_L999_R1_001.fastq.gz  
│	 └── ... (continued)  
├── illumiprocessor.conf 
└── illumiprocessor.log

```
I then ran a quick quality check on these reads. There are a number of other tools you can use to evaluate the effect of trimming on read counts, but the following script gives you a fairly accurate reading.

```
cd dmux-merged-clean/

for i in *;
do
    phyluce_assembly_get_fastq_lengths --input $i/split-adapter-quality-trimmed/ --csv;
done
```
Output:

```
All files in dir with aristochromis-christiae-ID25-READ2.fastq.gz,2778503,384269509,138.300915637,0.0142673237559,40,151,151.0
All files in dir with astatotilapia-burtoni-ID32-READ2.fastq.gz,2186967,298744995,136.602424728,0.0166710302187,40,151,151.0
All files in dir with aulonocara-stuartgranti-ID26-READ2.fastq.gz,2090484,288888592,138.192204293,0.0166037026004,40,151,151.0
All files in dir with bathybates-minor-ID90-READ2.fastq.gz,2213477,303924012,137.30615317,0.0163673516309,40,151,151.0
All files in dir with boulengerochromis-microlepis-ID79-READ2.fastq.gz,2312481,324622737,140.378553164,0.0147109324642,40,151,151.0
... (continued)
```
Each sample averaged about 2.5 million reads after illumiprocessing. Not bad. The results are stored in [fastqreads.csv] (fastqreads.csv).

###<a name="assembly">STEP 2: ASSEMBLE READS AND MATCH PROBES</a>
<a href="https://github.com/faircloth-lab/phyluce">phyluce</a> contains packages for three assembly methods: <a href="http://www.ebi.ac.uk/~zerbino/velvet/">velvet</a>, <a href="http://www.bcgsc.ca/platform/bioinfo/software/abyss">abyss</a>, and <a href="http://trinityrnaseq.sourceforge.net/">trinity</a>. The configuration files you feed into each assembly are similar, so you can make easy modifications to adapt to a different program. I use <a href="http://sourceforge.net/projects/trinityrnaseq/files/PREV_CONTENTS/previous_releases/">trinity</a> (vers. 2-25-2013) for this analysis. This older version generates significantly more contigs and runs much faster than the newer updates for complex reasons.

I created a config file that directed trinity to wherever my split-adapter-quality-trimmed reads are.

```
# make sure directory paths are correct. otherwise, trinity will not run.
aristochromis-christiae-ID25:/home/jimmyzheng/JZ_UCE_RUNS/dmux-merged-clean/aristochromis-christiae-ID25
astatotilapia-burtoni-ID32:/home/jimmyzheng/JZ_UCE_RUNS/dmux-merged-clean/astatotilapia-burtoni-ID32
aulonocara-stuartgranti-ID26:/home/jimmyzheng/JZ_UCE_RUNS/dmux-merged-clean/aulonocara-stuartgranti-ID26
bathybates-minor-ID90:/home/jimmyzheng/JZ_UCE_RUNS/dmux-merged-clean/bathybates-minor-ID90
... (continued)
```
I then saved this file in my `JZ_UCE_RUNS` base directory, from which the assembly will be generated. Trinity will create a new directory, which I specified below in the following script. The program runs at a program-defined k-mer value of 25. Smaller k-mer values allows for locating more overlapping sequences but smaller sized subsequences. Thus this can lead to more sequence ambiguities. Higher k-mer values leads to fewer ambiguities, longer contigs, but can take more computing time and space. You may adjust as you like, but from my experience (using a `min-kmer-coverage` of 60), the results were almost the same.
```
phyluce_assembly_assemblo_trinity \
    --config assembly.conf \
    --output /home/jimmyzheng/JZ_UCE_RUNS/trinity-assemblies \
    --clean \
    --cores 24
```
Once this finished (took about a day), I ran an additional quality check for the assemblies. The results are stored in [contigs.csv](contigs.csv).

```
for i in trinity-assemblies/contigs/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv;
done
```

Furthermore, it would be convenient to run a coverage check at this point, especially because these stats will be used in downstream analyses (particularly when [calling SNPs](#snps)). Before running this command, edit your `assembly.conf` by removing the "split-adapter-quality-trimmed" end, because you specify that in the command with `--subfolder`.

```
phyluce_assembly_get_trinity_coverage \
	--assemblies trinity-assemblies \
	--assemblo-config assembly.conf \
	--cores 12 \
	--clean \
	--subfolder split-adapter-quality-trimmed \
	--bwa-mem
```

Next up, we have to match our assembled contigs to our cichlid 1K UCE loci probe set. In the base directory, I web downloaded the file:

```
wget https://www.dropbox.com/s/0i7u8lf5cne2ljc/uce-fish-1k-probe-set.tar.gz?dl=0
```
A simple script will help us match the appropriate UCE loci sequences to our data. Be wary: most other model systems will require different probe sets - download the correct one, or else the number of UCE loci found will be pathetically low. You can find other probe sets from <a href="https://github.com/faircloth-lab/uce-probe-sets">Brant Faircloth's github</a>.

```
phyluce_assembly_match_contigs_to_probes \
    --contigs trinity-assemblies/contigs \
    --probes uce-1k-probes.fasta \
    --output cichlid-lastz
    # this is where all of the contigs.lastz files will be stored
```
.lastz files are aligned DNA files. What's most important in the `uce-search-results` directory is the probe.matches.sqlite database, which displays a table of all probe matches across species.

In order to access this probe.matches.sqlite database, I inputted the following code:

```
sqlite3 probe.matches.sqlite
sqlite> .output stdout  # display to the screen, not to file
sqlite> .headers on  # species names will be displayed
sqlite> .nullvalue --   # for non-matches, display a '--'
sqlite> SELECT * FROM matches LIMIT 10;  # limit to first 10 uce loci
```
The output is stored in [probe_matches.csv](uce-search-results/probe_matches.csv) but looks like the following:

```
# species names are abbreviated
uce         a_christiae		a_burtonigenus	a_stuartgranti
----------  --------------  --------------  --------------
uce-500     1               --              --
uce-501     1               --              --
uce-502     1               --              --
uce-503     1               1				1
uce-504     1               --              --
uce-505     1               --              --
uce-506     --              --              --
uce-507     1               --              --
uce-508     1               1               --
uce-509     1               1               1
... (continued)

```
###<a name="outgroup">STEP 3: INCLUDING OUTGROUP OR EXTERNAL DATA</a>
I also wanted to include **outgroup or external data** that is not contained within my assemblies. To do this, you have to assemble those raw data (if necessary), parse the assemblies (especially if not done by trinity or velvet), and then match the probes to the contigs. 

In my case, <a href="https://www.broadinstitute.org/ftp/pub/assemblies/fish/P_nyererei/PunNye1.0/">*pundamilia nyererei*</a>, as well as the outgroups <a href="http://ftp.ensembl.org/pub/release-82/fasta/oreochromis_niloticus/dna/">*oreochromis niloticus*</a> and <a href="https://www.broadinstitute.org/ftp/pub/assemblies/fish/N_brichardi/NeoBri1.0/">*neolamprologus brichardi*</a>, were publicly available assembled genomes, so I just needed to parse and match. Upon web downloading the data, I wrote the following to create an `outgroup-taxon-sets` directory that contains the uce-search-results for just the external data.

```
# to convert the fastas to 2bit files, processable by the next program
faToTwoBit H_burtoni_v1.assembly.fasta hapBur1.2bit
faToTwoBit N_brichardi_v1.assembly.fasta neoBri1.2bit
faToTwoBit Oreochromis_niloticus.Orenil1.0.dna.toplevel.fa oreNil1.2bit
faToTwoBit P_nyererei_v1.assembly.fasta punNye1.2bit

# separate directories for each external dataset
mkdir hapBur1; mkdir neoBri1; mkdir oreNil1; mkdir punNye1;

mkdir -p new-outgroups/outgroup-lastz; cd new-outgroups;

phyluce_probe_run_multiple_lastzs_sqlite \
	--db uce-1k-probes.sqlite \
	--output outgroup-lastz \
    --probefile uce-1k-probes.fasta \
    --scaffoldlist punNye1 oreNil1 hapBur1 neoBri1 \
    --genome-base-path . \
    --cores 12
```

This will generate DNA alignments for each external genome. In order to extract the fasta sequences from these lastz files for probe matching, we will need a config file to specify where the original external data are located.

```
#Config File: genome-sequence-location.conf in new-outgroups directory
[chromos]
oreochromis_niloticus_genome:/home/jimmyzheng/JZ_UCE_RUNS/outgroup-taxon-sets/oreochromis_niloticus_genome/oreochromis_niloticus_genome.2bit

[scaffolds]
p_nyererei:/home/jimmyzheng/JZ_UCE_RUNS/outgroup-taxon-sets/p_nyererei/p_nyererei.2bit
n_brichardi:/home/jimmyzheng/JZ_UCE_RUNS/outgroup-taxon-sets/n_brichardi/n_brichardi.2bit
```

```
phyluce_probe_slice_sequence_from_genomes \
	--conf ../genome-sequence-location.conf \
	--lastz outgroup-lastz \
	--output outgroup-fasta \
	--flank=1000 \
	--name-pattern "uce-1k-probes.fasta_v_{}.lastz.clean"
```
Once the fastas for each genome are sliced, I then copied these fastas over to the `trinity-assemblies/contigs` directory and matched contigs to probes, which generated the `cichlid-lastz` subdirectory, as described in the previous section.

From these lastz files, we now need to grab the actual contigs for further comparative analyses. This requires a config file called [taxon-set.conf] (taxon-set.conf), where I indicate what taxa goes into the final run.

```
# Config file: taxon-set.conf in base directory
[extended]
aristochromis-christiae-ID25
astatotilapia-burtoni-ID32
aulonocara-stuartgranti-ID26
capidochromis-eucinostomus-ID12
cheilotilapia-euchilus-ID23
cheilotilapia-rhodessi-ID24
copadichromis-trimaculatus-ID18
cyathichromis-obliquidens-ID3
docimodus-evelynae-ID10
fossorichromis-rostratus-ID28
genyochromis-mento-ID7
hemitilapia-oxyrhynchus-ID15
labeotropheus-fullbornei-ID8
labeotropheus-trewavase-ID6
labidochromis-gigas-ID30
melanochromis-auratus-ID4
melanochromis-kaskazini-ID33
metriaclima-greshakae-ID11
metriaclima-patricki-ID22
nimbochromis-polystigma-ID2
otopharynx-heterodon-ID13
otopharynx-lithobates-ID16
otopharynx-picta-ID14
petrotilapia-nigra-ID9
placidochromis-electra-ID21
placidochromis-milomo-ID17
psedotropheus-crabro-ID27
pseudotropheus-flavus-ID19
rhamphochromis-longiceps-ID81
simochromis-babaulti-ID48
stigmatochromis-woodei-ID20
taeniolethrinops-preorbitalis-ID29
tropheops-microstoma-ID1
tyrannochromis-nigriventer-ID5
hapBur1
neoBri1
oreNil1
punNye1
```
```
mkdir taxon-sets/final/
phyluce_assembly_get_match_counts \
	--locus-db cichlid-lastz/probe.matches.sqlite \
    --taxon-list-config taxon-set.conf \
    --taxon-group 'final' \
    --incomplete-matrix \
    --output taxon-sets/final/all-taxa-incomplete.conf 
```

The last step before alignments and trimming is to create an incomplete matrix for **all** of our data, including external genomes. I encountered some issues with naming, as dashes can choke the following program. The safest naming scheme, in my experience, is: `genus_species_###`, with underscores and numbers.

```
cd taxon-sets/final; mkdir log;
phyluce_assembly_get_fastas_from_match_counts \
    --contigs ../../trinity-assemblies/contigs \
    --locus-db ../../cichlid-lastz/probe.matches.sqlite \
    --match-count-output all-taxa-incomplete.conf \
    --output all-taxa-incomplete.fasta \
    --incomplete-matrix all-taxa-incomplete.incomplete \
    --log-path log
```
At this point, let's do another coverage check. We'll need to compare coverage across taxa if we choose to call SNPs. It's also just good to see how well your data represents this model system. To see a concise summary of coverage, input the following command:

```
phyluce_assembly_parse_trinity_coverage_log \
	--log ../../log/phyluce_assembly_get_trinity_coverage.log \
	--output malawi-contig-coverage.info.txt

# This output will be most important
phyluce_assembly_parse_trinity_coverage_for_uce_loci_log \
	--log ~/JZ_UCE_RUNS/taxon-sets/final/log/phyluce_assembly_get_trinity_coverage_for_uce_loci.log \
	--output malawi-UCE-coverage.info.txt
```

###<a name="alignment">STEP 4: ALIGNING SEQUENCES AND TRIMMING</a>

Once the incomplete data matrix is created, we need to explode the matrix into fasta files for each species based on each UCE loci sequence. To do this, I ran the following script:

```
cd taxon-sets/final

phyluce_assembly_explode_get_fastas_file \
    --input all-taxa-incomplete.fasta \
    --output-dir exploded-fastas \
    --by-taxon
    
for i in exploded-fastas/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv;
done

```
The summary file can be found in `taxon-sets/extended/` as [fasta_summary.csv] (fasta_summary.csv):

```
samples,contigs,total bp,mean length,95 CI length,min length,max length,median length,contigs >1kb
aristochromis-christiae-ID25.unaligned.fasta,1044,723419,692.930076628,4.21509495358,221,1362,700.0,12
aulonocara-stuartgranti-ID26.unaligned.fasta,1044,677615,649.05651341,5.03886400045,216,3831,652.0,4
bathybates-minor-ID90.unaligned.fasta,1045,698329,668.257416268,4.30500853163,229,1017,683.0,2
capidochromis-eucinostomus-ID12.unaligned.fasta,1027,776873,756.448880234,5.30036625092,201,2990,772.0,45
cheilotilapia-euchilus-ID23.unaligned.fasta,1042,718002,689.061420345,4.24679134664,246,1231,696.0,8
... (continued)
```
Now for **aligning the data**. Be sure to specify the right number of taxa, *including external and outgroup taxa*. Otherwise, you might end up creating an extremely incomplete matrix, affecting downstream analyses. Send all log files into the `log` subdirectory.

```
phyluce_align_seqcap_align \
    --fasta all-taxa-incomplete.fasta \
    --output mafft-nexus-edge-trimmed \
    --taxa 38 \
    --aligner mafft \
    --cores 12 \
    --incomplete-matrix \
    --log-path log
```
What the seqcap align program does is both align all sequences and trim ragged ends if present. When we create the alignment, there may be a bunch of bases at the end that hang over. What we are doing here with edge trimming is pruning and cleaning up our data.

Time to see what our alignment looks like.

```
phyluce_align_get_align_summary_data \
    --alignments mafft-nexus-edge-trimmed \
    --cores 12 \
    --log-path log
```
The log file can be found in [phyluce_align_get_align_summary_data.log] (taxon-sets/extended/log/phyluce_align_get_align_summary_data.log), but I have also enclosed it here:

```
2015-12-07 20:46:06,866 - phyluce_align_get_align_summary_data - INFO - ========= Starting phyluce_align_get_align_summary_data =========
2015-12-07 20:46:06,866 - phyluce_align_get_align_summary_data - INFO - Version: 1.5.0
2015-12-07 20:46:06,866 - phyluce_align_get_align_summary_data - INFO - Argument --alignments: /home/jimmyzheng/JZ_UCE_RUNS/taxon-sets/final/mafft-untrimmed-gblocks
2015-12-07 20:46:06,866 - phyluce_align_get_align_summary_data - INFO - Argument --cores: 12
2015-12-07 20:46:06,867 - phyluce_align_get_align_summary_data - INFO - Argument --input_format: nexus
2015-12-07 20:46:06,867 - phyluce_align_get_align_summary_data - INFO - Argument --log_path: /home/jimmyzheng/JZ_UCE_RUNS/taxon-sets/final/log
2015-12-07 20:46:06,867 - phyluce_align_get_align_summary_data - INFO - Argument --show_taxon_counts: False
2015-12-07 20:46:06,867 - phyluce_align_get_align_summary_data - INFO - Argument --verbosity: INFO
2015-12-07 20:46:06,867 - phyluce_align_get_align_summary_data - INFO - Getting alignment files
2015-12-07 20:46:06,874 - phyluce_align_get_align_summary_data - INFO - Computing summary statistics using 12 cores
2015-12-07 20:46:09,070 - phyluce_align_get_align_summary_data - INFO - ----------------------- Alignment summary -----------------------
2015-12-07 20:46:09,070 - phyluce_align_get_align_summary_data - INFO - [Alignments] loci:	1,086
2015-12-07 20:46:09,070 - phyluce_align_get_align_summary_data - INFO - [Alignments] length:	641,561
2015-12-07 20:46:09,070 - phyluce_align_get_align_summary_data - INFO - [Alignments] mean:	590.76
2015-12-07 20:46:09,070 - phyluce_align_get_align_summary_data - INFO - [Alignments] 95% CI:	8.46
2015-12-07 20:46:09,071 - phyluce_align_get_align_summary_data - INFO - [Alignments] min:	172
2015-12-07 20:46:09,071 - phyluce_align_get_align_summary_data - INFO - [Alignments] max:	2,169
2015-12-07 20:46:09,072 - phyluce_align_get_align_summary_data - INFO - ------------------------- Taxon summary -------------------------
2015-12-07 20:46:09,072 - phyluce_align_get_align_summary_data - INFO - [Taxa] mean:		36.28
2015-12-07 20:46:09,072 - phyluce_align_get_align_summary_data - INFO - [Taxa] 95% CI:	0.31
2015-12-07 20:46:09,072 - phyluce_align_get_align_summary_data - INFO - [Taxa] min:		3
2015-12-07 20:46:09,072 - phyluce_align_get_align_summary_data - INFO - [Taxa] max:		38
2015-12-07 20:46:09,073 - phyluce_align_get_align_summary_data - INFO - ----------------- Missing data from trim summary ----------------
2015-12-07 20:46:09,073 - phyluce_align_get_align_summary_data - INFO - [Missing] mean:	0.00
2015-12-07 20:46:09,073 - phyluce_align_get_align_summary_data - INFO - [Missing] 95% CI:	0.00
2015-12-07 20:46:09,073 - phyluce_align_get_align_summary_data - INFO - [Missing] min:	0.00
2015-12-07 20:46:09,073 - phyluce_align_get_align_summary_data - INFO - [Missing] max:	0.00
2015-12-07 20:46:09,085 - phyluce_align_get_align_summary_data - INFO - -------------------- Character count summary --------------------
2015-12-07 20:46:09,085 - phyluce_align_get_align_summary_data - INFO - [All characters]	23,119,842
2015-12-07 20:46:09,085 - phyluce_align_get_align_summary_data - INFO - [Nucleotides]		22,748,282
2015-12-07 20:46:09,087 - phyluce_align_get_align_summary_data - INFO - ---------------- Data matrix completeness summary ---------------
2015-12-07 20:46:09,087 - phyluce_align_get_align_summary_data - INFO - [Matrix 50%]		1059 alignments
2015-12-07 20:46:09,087 - phyluce_align_get_align_summary_data - INFO - [Matrix 55%]		1054 alignments
2015-12-07 20:46:09,087 - phyluce_align_get_align_summary_data - INFO - [Matrix 60%]		1052 alignments
2015-12-07 20:46:09,088 - phyluce_align_get_align_summary_data - INFO - [Matrix 65%]		1049 alignments
2015-12-07 20:46:09,088 - phyluce_align_get_align_summary_data - INFO - [Matrix 70%]		1046 alignments
2015-12-07 20:46:09,088 - phyluce_align_get_align_summary_data - INFO - [Matrix 75%]		1041 alignments
2015-12-07 20:46:09,088 - phyluce_align_get_align_summary_data - INFO - [Matrix 80%]		1026 alignments
2015-12-07 20:46:09,088 - phyluce_align_get_align_summary_data - INFO - [Matrix 85%]		1020 alignments
2015-12-07 20:46:09,088 - phyluce_align_get_align_summary_data - INFO - [Matrix 90%]		990 alignments
2015-12-07 20:46:09,088 - phyluce_align_get_align_summary_data - INFO - [Matrix 95%]		933 alignments
2015-12-07 20:46:09,088 - phyluce_align_get_align_summary_data - INFO - ------------------------ Character counts -----------------------
2015-12-07 20:46:09,089 - phyluce_align_get_align_summary_data - INFO - [Characters] '-' is present 371,560 times
2015-12-07 20:46:09,089 - phyluce_align_get_align_summary_data - INFO - [Characters] 'A' is present 6,345,309 times
2015-12-07 20:46:09,089 - phyluce_align_get_align_summary_data - INFO - [Characters] 'C' is present 4,997,282 times
2015-12-07 20:46:09,089 - phyluce_align_get_align_summary_data - INFO - [Characters] 'G' is present 5,045,547 times
2015-12-07 20:46:09,089 - phyluce_align_get_align_summary_data - INFO - [Characters] 'T' is present 6,360,144 times
2015-12-07 20:46:09,089 - phyluce_align_get_align_summary_data - INFO - ========= Completed phyluce_align_get_align_summary_data ========
```
Our 75% matrix looks good at 1041 alignments and our 95% matrix has a whopping 933. We can also choose not to trim and see how our alignments look. In our case, they do not look much different, so we move forward. The script to run with just internal trimming is as follows:

```
phyluce_align_seqcap_align \
    --fasta all-taxa-incomplete.fasta \
    --output mafft-nexus-internal-trimmed \
    --taxa 38 \
    --aligner mafft \
    --cores 12 \
    --incomplete-matrix \
    --output-format fasta \
    --no-trim \ # MUST SPECIFY to prevent trimming
    --log-path log
```
Then pass alignments to Gblocks wrapper and remove UCE loci names from the nexus.

```
phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
    --alignments mafft-nexus-internal-trimmed \
    --output mafft-nexus-internal-trimmed-gblocks \
    --b1 0.5 \
    --b4 8 \
    --cores 12 \
    --log log

phyluce_align_remove_locus_name_from_nexus_lines \
    --alignments mafft-nexus-internal-trimmed-gblocks \
    --output mafft-nexus-internal-trimmed-gblocks-clean \
    --cores 12 \
    --log-path log
```
You can also directly remove locus names from nexus lines from the edge-trimmed alignments. I chose to go with the internally trimmed alignment, though there shouldn't be any significant difference.

##<a name="phylogeny">Final Matrices and  Phylogenetic Analysis</a>
Once I finished the data preparation, I moved on to generating the final matrix. I chose to go with the 75% incomplete matrix for RAxML and ExaBayes analyses. To do that, I ran:

```
phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-nexus-internal-trimmed-gblocks-clean \
    --taxa 38 \ # make sure this number has been accurate throughout
    --percent 0.95 \ # to specify 75% incompleteness
    --output mafft-nexus-internal-trimmed-gblocks-clean-75p \
    --cores 12 \
    --log-path log

# FOR RAXML    
phyluce_align_format_nexus_files_for_raxml \
    --alignments mafft-nexus-internal-trimmed-gblocks-clean-75p \
    --output mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml \
    --charsets \
    --log-path log
```

My directory structure, shortened, looks like this:

```
JZ_UCE_RUNS
├── taxon-sets
│	 └── extended 
│	     ├── all-taxa-incomplete.conf  
│	     ├── all-taxa-incomplete.fasta  
│	     ├── all-taxa-incomplete.incomplete  
│	     ├── exploded-fastas  
│	     ├── mafft-nexus-edge-trimmed 
│	     ├── mafft-nexus-internal-trimmed  
│	     ├── mafft-nexus-internal-trimmed-gblocks  
│	     ├── mafft-nexus-internal-trimmed-gblocks-clean  
│	     ├── mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml-part  
│	     ├── mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml-unpart  
│	     ├── mafft-nexus-internal-trimmed-gblocks-clean-75p-exabayes-part  
│	     ├── mafft-nexus-internal-trimmed-gblocks-clean-75p-exabayes-unpart  
│	     └── log
└── ... (continued)
```

Then, I began my RAxML analyses. 

```
raxmlHPC-PTHREADS-SSE3 \
    -m GTRGAMMA \
    -N 20 \
    -p 18262063 \
    -n best \
    -s mafft-gblocks-95p.phylip \
    -T 12 \
    -o oreNil1
```
Compute bootstraps:

```
raxmlHPC-PTHREADS-SSE3 \
    -m GTRGAMMA \
    -N autoMRE \
    -p 172421953 \
    -b 92815632 \
    -n bootreps \
    -s mafft-gblocks-95p.phylip \
    -o oreNil1 \
    -T 12
```
Reconcile tree and support values:

```
raxmlHPC-SSE3 \  
    -m GTRGAMMA \  
    -n bestML.bootstrap\
    -f b \
    -t RAxML_bestTree.best \
    -z RAxML_bootstrap.bootreps
```
Trees are stored in `JZ-UCE-RUNS/trees/new-assembly`. I am currently generating the ExaBayes trees for both 75p and 95p.

The commands are as follows:
```
mpirun -np 16 exabayes -f mafft-gblocks-95p.phylip -s 113012905 -n run1 -c config.nexus -q aln.part -R 4 -C 2  

consense -f ExaBayes_topologies.run1.* -n malawis
postProcParam -f ExaBayes_parameters.run1.* -n malawis
```

To run the partitioned analyses, you need to use PartitionFinder v1.1.0 by Robert Lanfear, using the `--rcluster` option and `GTR+GAMMA` model selection. The search for ideal partitions will take up to 3 days. Once those are done, input them into your RAxML and ExaBayes runs (`-q`).

To load in initial partitions to PartitionFinder, I nano'ed the file `mafft-nexus-internal-trimmed-gblocks-clean-75p.charsets` in the `mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml-part` directory.

```
begin sets;
charset 'uce-49.nexus' = 226277-226956;
charset 'uce-399.nexus' = 126906-127537;
charset 'uce-1160.nexus' = 370136-370624;
charset 'uce-453.nexus' = 328938-329516;
charset 'uce-1216.nexus' = 319363-319970;
charset 'uce-615.nexus' = 78810-79379;
charset 'uce-953.nexus' = 266557-267239;
charset 'uce-604.nexus' = 490391-490965;
... (continued)
```

PF will figure out which of these charsets are more similar to others and group them into larger partition subsets.

Currently running the partitions...

##<a name="snps">Calling SNPs and Generating Species Trees</a>

For extremely fast-changing model systems, it might be a good idea to call SNPs to zoom in on the few but important polymorphisms that do affect speciation analyses. In order to do this, we use a different method altogether.

We'll need to grab one sample that has the highest coverage. We can look into our `malawi-UCE-coverage.info.txt` file and sort by coverage. You can do this in Excel.

Then we get match counts and fastas for JUST that one species - in my case, it was tyrannochromis-nigriventer-ID5.

```
phyluce_assembly_get_match_counts \
    --locus-db ../../cichlid-lastz/probe.matches.sqlite \
    --taxon-list-config ../../taxon-set.conf \
    --taxon-group 'one' \
    --output tyrannochromis-nigriventer.conf
        
phyluce_assembly_get_fastas_from_match_counts \
    --contigs ../../trinity-assemblies/contigs \
    --locus-db ../../cichlid-lastz/probe.matches.sqlite \
    --match-count-output tyrannochromis-nigriventer.conf \
    --output tyrannochromis-nigriventer.fasta

# index into .bam files
bwa index tyrannochromis-nigriventer.fasta
mkdir reference/
mv * reference/
```
Then, you're going to create a config file for the auto BWA-runner.

```
# Config File: tyr-nig-snps.conf

[reference]
tyrannochromis-nigriventer.fasta:/home/jimmyzheng/JZ_UCE_RUNS/taxon-sets/final/references

[individuals]
aristochromis-christiae-ID25:/home/jimmyzheng/JZ_UCE_RUNS/dmux-merged-clean/aristochromis-christiae-ID25
astatotilapia-burtoni-ID32:/home/jimmyzheng/JZ_UCE_RUNS/dmux-merged-clean/astatotilapia-burtoni-ID32
aulonocara-stuartgranti-ID26:/home/jimmyzheng/JZ_UCE_RUNS/dmux-merged-clean/aulonocara-stuartgranti-ID26
bathybates-minor-ID90:/home/jimmyzheng/JZ_UCE_RUNS/dmux-merged-clean/bathybates-minor-ID90
...

[flowcell]
aristochromis-christiae-ID25:H13T7BGXX
astatotilapia-burtoni-ID32:H13T7BGXX
aulonocara-stuartgranti-ID26:H13T7BGXX
bathybates-minor-ID90:H13T7BGXX
...
```
```
phyluce_snp_bwa_align \
	--config tyr-nig-snps.conf \
	--output ../ \
	--cores 12 \
	--subfolder split-adapter-quality-trimmed \
	--mem
```
This essentially creates BAM files for each taxon. Next, we have to align and merge them into one `cichlid-merged.bam` file alignment.

```
java -Xmx40g -jar ~/anaconda/jar/MergeSamFiles.jar \
	I=../aristochromis-christiae-ID25/aristochromis-christiae-ID25-CL-RG-MD-M.bam \
	I=../astatotilapia-burtoni-ID32/astatotilapia-burtoni-ID32-CL-RG-MD-M.bam \
	I=../aulonocara-stuartgranti-ID26/aulonocara-stuartgranti-ID26-CL-RG-MD-M.bam \
	I=../bathybates-minor-ID90/bathybates-minor-ID90-CL-RG-MD-M.bam \
	... (continued)
	SO=coordinate \
    AS=true \
    VALIDATION_STRINGENCY=LENIENT \
    O=cichlids-merged.bam

samtools index cichlids-merged.bam #index the BAM file
```
You'll want to create a sequence dictionary for your one reference taxon.

```
java -Xmx2g -jar ~/anaconda/jar/CreateSequenceDictionary.jar \
R=tyrannochromis-nigriventer.fasta \
O=tyrannochromis-nigriventer.dict

samtools faidx tyrannochromis-nigriventer.fasta
```
Now, you'll find the indels, which you will filter out SNPs from later.

```
java -Xmx20g -jar ~/anaconda/jar/GenomeAnalysisTK.jar \
	-T RealignerTargetCreator \
    -R tyrannochromis-nigriventer.fasta \
    -I cichlids-merged.bam \
    --minReadsAtLocus 4 \
    -o cichlids-merged.intervals \
    -nt 12 
```

Realign the BAM based on these indel intervals:

```
java -Xmx20g -jar ~/anaconda/jar/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R tyrannochromis-nigriventer.fasta \
    -I cichlids-merged.bam \
    -targetIntervals cichlids-merged.intervals \
    -LOD 3.0 \
    -o cichlids-realigned.bam
```
Call the SNPs:

```
java -Xmx20g -jar ~/anaconda/jar/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R tyrannochromis-nigriventer.fasta \
    -I cichlids-realigned.bam \
    -gt_mode DISCOVERY \
    -stand_call_conf 30 \
    -stand_emit_conf 10 \
    -o cichlids-rawSNPS-Q30.vcf
```
Call the indels:

```
java -Xmx20g -jar ~/anaconda/jar/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R tyrannochromis-nigriventer.fasta \
    --max_alternate_alleles 15 \ 
    #default is 6, value should depend on analysis
    -I cichlids-realigned.bam \
    -gt_mode DISCOVERY \
    -glm INDEL \
    -stand_call_conf 30 \
    -stand_emit_conf 10 \
    -o cichlids-inDels-Q30.vcf
```
Filter based on poor quality SNPs (within 5 bp of an indel, SNPs in clusters > size 10, SNP loci with quality < 30, SNP loci with QD values below 2 and SNP loci failing a filtering formula `MQ0 >= 4 && ((MQ0/(1.0 * DP)) > 0.1)`).

```
java -Xmx20g -jar ~/anaconda/jar/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R tyrannochromis-nigriventer.fasta \
    -V cichlids-rawSNPS-Q30.vcf \
    --mask cichlids-inDels-Q30.vcf \
    --maskExtension 5 \
    --maskName InDel \
    --clusterWindowSize 10 \
    --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
    --filterName "Bad Validation" \
    --filterExpression "QUAL < 25.0" \
    --filterName "LowQual" \
    --filterExpression "QD < 2.0" \
    --filterName "Low Variant Confidence" \
    -o cichlids-Q30-QD2-LOW-STRICT.vcf
```

Remove these poorly validated SNPs using vcftools. This will output the number of SNPs removed and kept. The output file will be used to create the actual species tree.

```
/home/jimmyzheng/vcftools/src/cpp/vcftools \
	--vcf cichlids-Q30-QD2-LOW-STRICT.vcf \
	--remove-filtered-all \
	--max-missing 0.75 \
	--recode \
	--out Q30-QD2-MISS_0.75
```
You'll want to create a nexus file to read into BEAST v2.2.1 and SNAPP, one of its packages. SNAPP will create the species tree.

```
# This text file can give you a look at the actual alignment
phyluce_snp_convert_vcf_to_structure \
	--input Q30-QD2-MISS_0.75.recode.vcf \
	--output Q30-QD2-MISS_0.75.recode.structure.txt \
	--filter-informative

phyluce_snp_convert_vcf_to_snapp \
	--input Q30-QD2-MISS_0.75.recode.vcf \
	--output Q30-QD2-MISS_0.75.recode.snapp.nexus \
	--filter-informative

# Append under begin nexus:
# datatype=standard symbols="012"
```
In order to run SNAPP, you have to convert your binary nexus file into a .xml file readable by BEAST. Use BEAUti as part of the v2.2.1 patch to generate that. Calculate forward and backward mutation rates and leave the priors as is. Then run into BEAST using this command:

```
~/jdk1.8.0_65/bin/java -jar ~/beast/lib/beast.jar -threads 24 malawis.xml
```
Note: you must use Java 1.8+ versions in order for BEAST v2 to work. I downloaded the patch locally.

Currently running the species tree analysis...