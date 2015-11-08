#Ultraconserved Element (UCE) Analyses for Malawi Cichlids
This github serves the purpose of running analyses on DNA ultra-conserved elements data for Malawi cichlids. You will find configuration files, pipeline scripts, and reconstructed trees.
Credit goes to Brant Faircloth for constructing the <a href="https://github.com/faircloth-lab/phyluce">phyluce pipeline</a>.


***Jimmy Zheng***  
*November 7, 2015*  

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
Next, I created a configuration file (.conf) for the illumiprocessor program, which will copy the raw reads into a different directory and remove contaminating adapter sequences and low-quality bases. This is where I match up each species to its respective adapter identifiers. See below. The following file includes both Tanganyika and Malawi cichlids.

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
i7-07_06:GTCACTGT
i7-07_07:CCGTAAGA
i7-07_08:GTCTTGCA
i7-07_09:ACTGAGGT
i7-07_10:CACTGACA
i7-07_11:CAGTCCAA
i7-07_12:TCGACATC
i7-08_01:GAGTCTCT
i7-08_02:ATTGCGTG
i7-08_03:TGCGAACT
i7-08_04:GGTAGTGT
i7-08_05:CAACGGAT
i7-08_06:TACGGTTG
i7-08_07:CAAGTGCA
i7-08_08:ATGCACGA
i7-08_09:AGCAAGCA
i7-08_10:CTAGGTGA
i7-08_11:ATCGCCAT
i7-08_12:TCGAACCA
i7-09_05:GATAGCGA
i7-09_06:CCAAGACT
i7-09_07:GGATACCA
i7-09_08:CGTAGGTT
i7-09_09:CTGGTTCT
i7-09_10:ATGCCAAC
i7-09_11:ATATGCGC
i5-02_H:GCTTCGAA
i5-03_A:GTGGTGTT
i5-03_B:ACAGCTCA
i5-03_C:TTCCTGTG
i5-03_D:GGTTGTCA
i5-03_E:ACGGAACA
i5-03_F:TCTCTAGG
i5-03_G:CGTTATGC
i5-03_H:AAGCACTG
i5-04_A:GAGATACG

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
cheilotilapia-rhodessi-ID24_S153:i7-08_05,i5-03_H
copadichromis-trimaculatus-ID18_S169:i7-08_06,i5-03_A
cyathichromis-obliquidens-ID3_S156:i7-07_03,i5-02_H
cyprichromis-pavo-ID84_S237:i7-09_11,i5-04_A
docimodus-evelynae-ID10_S162:i7-07_10,i5-02_H
ectodus-descampsii-ID55_S200:i7-08_08,i5-03_D
eretmodus-cyanostictus-ID82_S222:i7-08_11,i5-03_F
fossorichromis-rostratus-ID28_S173:i7-08_04,i5-03_B
genyochromis-mento-ID7_S159:i7-07_07,i5-02_H
gnathochromis-permaxillaris-ID60_S230:i7-08_11,i5-03_H
haplochromis-thereuterion-ID71_S214:i7-08_12,i5-03_E
haplotaxodon-microlepis-ID77_S218:i7-08_06,i5-03_F
harpagochromis-golden-duck-ID73_S216:i7-08_02,i5-03_F
harpagochromis-orange-rock-hunter-ID89_S227:i7-08_06,i5-03_G
harpagochromis-serranus-ID87_S231:i7-09_05,i5-04_A
harpagochromis-two-stripe-white-lip-ID88_S226:i7-08_05,i5-03_G
hemitilapia-oxyrhynchus-ID15_S167:i7-08_03,i5-03_A
julidochromis-dickfeldi-ID45_S189:i7-08_09,i5-03_C
labeotropheus-fullbornei-ID8_S160:i7-07_08,i5-02_H
labeotropheus-trewavase-ID6_S158:i7-07_06,i5-02_H
labidochromis-gigas-ID30_S175:i7-08_06,i5-03_B
lamprologus-callipterus-ID56_S201:i7-08_09,i5-03_D
lamprologus-signatus-ID70_S213:i7-08_11,i5-03_E
lamprologus-speciosus-ID65_S209:i7-08_06,i5-03_E
lepidiolamprologus-attenuatus-ID47_S191:i7-08_11,i5-03_C
lepidiolamprologus-boulengeri-ID36_S181:i7-08_12,i5-03_B
lepidiolamprologus-hecqui-ID58_S203:i7-08_11,i5-03_D
lepidiolamprologus-leimairii-ID69_S233:i7-09_07,i5-04_A
lepidiolamprologus-meeli-kipili-ID40_S229:i7-08_08,i5-03_G
lepidiolamprologus-nkambae-ID41_S185:i7-08_05,i5-03_C
lepidiolamprologus-profundicola-ID76_S235:i7-09_09,i5-04_A
lepidolamprologus-elongatus-ID34_S179:i7-08_10,i5-03_B
limnochromis-auritus-ID38_S183:i7-08_02,i5-03_C
lipochromis-parvidens-ID74_S234:i7-09_08,i5-04_A
melanochromis-auratus-ID4_S157:i7-07_04,i5-02_H
melanochromis-kaskazini-ID33_S178:i7-08_09,i5-03_B
metriaclima-greshakae-ID11_S163:i7-07_11,i5-02_H
metriaclima-patricki-ID22_S151:i7-08_03,i5-03_H
neolamprologus-bifasciatus-ID80_S236:i7-09_10,i5-04_A
neolamprologus-brichardi-ID50_S194:i7-08_02,i5-03_D
neolamprologus-cylindricus-ID68_S212:i7-08_09,i5-03_E
neolamprologus-fasciatus-ID86_S225:i7-08_03,i5-03_G
neolamprologus-furcifer-ID66_S210:i7-08_07,i5-03_E
neolamprologus-helianthus-ID83_S223:i7-08_12,i5-03_F
neolamprologus-longicaudata-ID54_S199:i7-08_07,i5-03_D
neolamprologus-longior-ID35_S180:i7-08_11,i5-03_B
neolamprologus-modestus-ID61_S205:i7-08_02,i5-03_E
neolamprologus-niger-ID59_S204:i7-08_12,i5-03_D
neolamprologus-nigriventer-ID67_S211:i7-08_08,i5-03_E
neolamprologus-obscurus-ID53_S198:i7-08_06,i5-03_D
neolamprologus-olivaceous-ID85_S224:i7-08_02,i5-03_G
neolamprologus-prochilus-ID51_S195:i7-08_03,i5-03_D
neolamprologus-pulcher-ID42_S186:i7-08_06,i5-03_C
neolamprologus-savoryi-ID39_S184:i7-08_03,i5-03_C
neolamprologus-sexfasciatus-ID63_S207:i7-08_04,i5-03_E
neolamprologus-tetracanthus-ID46_S190:i7-08_10,i5-03_C
neolamprologus-walteri-ID52_S196:i7-08_04,i5-03_D
nimbochromis-polystigma-ID2_S155:i7-07_02,i5-02_H
ophthalmotilapia-boops-ID44_S188:i7-08_08,i5-03_C
ophthalmotilapia-nasuta-ID37_S182:i7-08_01,i5-03_C
otopharynx-heterodon-ID13_S165:i7-08_01,i5-03_A
otopharynx-lithobates-ID16_S168:i7-08_04,i5-03_A
otopharynx-picta-ID14_S166:i7-08_02,i5-03_A
petrotilapia-nigra-ID9_S161:i7-07_09,i5-02_H
placidochromis-electra-ID21_S197:i7-08_05,i5-03_D
placidochromis-milomo-ID17_S232:i7-09_06,i5-04_A
prognathochromis-perrieri-ID49_S193:i7-08_01,i5-03_D
psedotropheus-crabro-ID27_S172:i7-08_03,i5-03_B
pseudotropheus-flavus-ID19_S149:i7-08_01,i5-03_H
ptyochromis-cf-ID72_S215:i7-08_01,i5-03_F
pundamila-pundamila-ID75_S217:i7-08_04,i5-03_F
pyxichromis-orthostoma-ID62_S206:i7-08_03,i5-03_E
rhamphochromis-longiceps-ID81_S221:i7-08_10,i5-03_F
simochromis-babaulti-ID48_S192:i7-08_12,i5-03_C
stigmatochromis-woodei-ID20_S150:i7-08_02,i5-03_H
taeniolethrinops-preorbitalis-ID29_S174:i7-08_05,i5-03_B
telmatochromis-dhonti-ID43_S187:i7-08_07,i5-03_C
telmatochromis-temporalis-ID64_S208:i7-08_05,i5-03_E
trematochromis-benthicola-ID78_S219:i7-08_07,i5-03_F
tropheops-microstoma-ID1_S154:i7-07_01,i5-02_H
tyrannochromis-nigriventer-ID5_S148:i7-06_01,i5-02_H
variabilichromis-moorii-ID31_S176:i7-08_07,i5-03_B

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
cheilotilapia-rhodessi-ID24_S153:cheilotilapia-rhodessi-ID24
copadichromis-trimaculatus-ID18_S169:copadichromis-trimaculatus-ID18
cyathichromis-obliquidens-ID3_S156:cyathichromis-obliquidens-ID3
cyprichromis-pavo-ID84_S237:cyprichromis-pavo-ID84
docimodus-evelynae-ID10_S162:docimodus-evelynae-ID10
ectodus-descampsii-ID55_S200:ectodus-descampsii-ID55
eretmodus-cyanostictus-ID82_S222:eretmodus-cyanostictus-ID82
fossorichromis-rostratus-ID28_S173:fossorichromis-rostratus-ID28
genyochromis-mento-ID7_S159:genyochromis-mento-ID7
gnathochromis-permaxillaris-ID60_S230:gnathochromis-permaxillaris-ID60
haplochromis-thereuterion-ID71_S214:haplochromis-thereuterion-ID71
haplotaxodon-microlepis-ID77_S218:haplotaxodon-microlepis-ID77
harpagochromis-golden-duck-ID73_S216:harpagochromis-golden-duck-ID73
harpagochromis-orange-rock-hunter-ID89_S227:harpagochromis-orange-rock-hunter-ID89
harpagochromis-serranus-ID87_S231:harpagochromis-serranus-ID87
harpagochromis-two-stripe-white-lip-ID88_S226:harpagochromis-two-stripe-white-lip-ID88
hemitilapia-oxyrhynchus-ID15_S167:hemitilapia-oxyrhynchus-ID15
julidochromis-dickfeldi-ID45_S189:julidochromis-dickfeldi-ID45
labeotropheus-fullbornei-ID8_S160:labeotropheus-fullbornei-ID8
labeotropheus-trewavase-ID6_S158:labeotropheus-trewavase-ID6
labidochromis-gigas-ID30_S175:labidochromis-gigas-ID30
lamprologus-callipterus-ID56_S201:lamprologus-callipterus-ID56
lamprologus-signatus-ID70_S213:lamprologus-signatus-ID70
lamprologus-speciosus-ID65_S209:lamprologus-speciosus-ID65
lepidiolamprologus-attenuatus-ID47_S191:lepidiolamprologus-attenuatus-ID47
lepidiolamprologus-boulengeri-ID36_S181:lepidiolamprologus-boulengeri-ID36
lepidiolamprologus-hecqui-ID58_S203:lepidiolamprologus-hecqui-ID58
lepidiolamprologus-leimairii-ID69_S233:lepidiolamprologus-leimairii-ID69
lepidiolamprologus-meeli-kipili-ID40_S229:lepidiolamprologus-meeli-kipili-ID40
lepidiolamprologus-nkambae-ID41_S185:lepidiolamprologus-nkambae-ID41
lepidiolamprologus-profundicola-ID76_S235:lepidiolamprologus-profundicola-ID76
lepidolamprologus-elongatus-ID34_S179:lepidolamprologus-elongatus-ID34
limnochromis-auritus-ID38_S183:limnochromis-auritus-ID38
lipochromis-parvidens-ID74_S234:lipochromis-parvidens-ID74
melanochromis-auratus-ID4_S157:melanochromis-auratus-ID4
melanochromis-kaskazini-ID33_S178:melanochromis-kaskazini-ID33
metriaclima-greshakae-ID11_S163:metriaclima-greshakae-ID11
metriaclima-patricki-ID22_S151:metriaclima-patricki-ID22
neolamprologus-bifasciatus-ID80_S236:neolamprologus-bifasciatus-ID80
neolamprologus-brichardi-ID50_S194:neolamprologus-brichardi-ID50
neolamprologus-cylindricus-ID68_S212:neolamprologus-cylindricus-ID68
neolamprologus-fasciatus-ID86_S225:neolamprologus-fasciatus-ID86
neolamprologus-furcifer-ID66_S210:neolamprologus-furcifer-ID66
neolamprologus-helianthus-ID83_S223:neolamprologus-helianthus-ID83
neolamprologus-longicaudata-ID54_S199:neolamprologus-longicaudata-ID54
neolamprologus-longior-ID35_S180:neolamprologus-longior-ID35
neolamprologus-modestus-ID61_S205:neolamprologus-modestus-ID61
neolamprologus-niger-ID59_S204:neolamprologus-niger-ID59
neolamprologus-nigriventer-ID67_S211:neolamprologus-nigriventer-ID67
neolamprologus-obscurus-ID53_S198:neolamprologus-obscurus-ID53
neolamprologus-olivaceous-ID85_S224:neolamprologus-olivaceous-ID85
neolamprologus-prochilus-ID51_S195:neolamprologus-prochilus-ID51
neolamprologus-pulcher-ID42_S186:neolamprologus-pulcher-ID42
neolamprologus-savoryi-ID39_S184:neolamprologus-savoryi-ID39
neolamprologus-sexfasciatus-ID63_S207:neolamprologus-sexfasciatus-ID63
neolamprologus-tetracanthus-ID46_S190:neolamprologus-tetracanthus-ID46
neolamprologus-walteri-ID52_S196:neolamprologus-walteri-ID52
nimbochromis-polystigma-ID2_S155:nimbochromis-polystigma-ID2
ophthalmotilapia-boops-ID44_S188:ophthalmotilapia-boops-ID44
ophthalmotilapia-nasuta-ID37_S182:ophthalmotilapia-nasuta-ID37
otopharynx-heterodon-ID13_S165:otopharynx-heterodon-ID13
otopharynx-lithobates-ID16_S168:otopharynx-lithobates-ID16
otopharynx-picta-ID14_S166:otopharynx-picta-ID14
petrotilapia-nigra-ID9_S161:petrotilapia-nigra-ID9
placidochromis-electra-ID21_S197:placidochromis-electra-ID21
placidochromis-milomo-ID17_S232:placidochromis-milomo-ID17
prognathochromis-perrieri-ID49_S193:prognathochromis-perrieri-ID49
psedotropheus-crabro-ID27_S172:psedotropheus-crabro-ID27
pseudotropheus-flavus-ID19_S149:pseudotropheus-flavus-ID19
ptyochromis-cf-ID72_S215:ptyochromis-cf-ID72
pundamila-pundamila-ID75_S217:pundamila-pundamila-ID75
pyxichromis-orthostoma-ID62_S206:pyxichromis-orthostoma-ID62
rhamphochromis-longiceps-ID81_S221:rhamphochromis-longiceps-ID81
simochromis-babaulti-ID48_S192:simochromis-babaulti-ID48
stigmatochromis-woodei-ID20_S150:stigmatochromis-woodei-ID20
taeniolethrinops-preorbitalis-ID29_S174:taeniolethrinops-preorbitalis-ID29
telmatochromis-dhonti-ID43_S187:telmatochromis-dhonti-ID43
telmatochromis-temporalis-ID64_S208:telmatochromis-temporalis-ID64
trematochromis-benthicola-ID78_S219:trematochromis-benthicola-ID78
tropheops-microstoma-ID1_S154:tropheops-microstoma-ID1
tyrannochromis-nigriventer-ID5_S148:tyrannochromis-nigriventer-ID5
variabilichromis-moorii-ID31_S176:variabilichromis-moorii-ID31

```

Then I ran <a href="http://illumiprocessor.readthedocs.org/en/latest/usage.html">illumiprocessor</a> using the following code:

```
illumiprocessor \
    --input raw-fastq/ \
    --output clean-fastq \
    --config illumiprocessor.conf \
    --cores 12
```


Your directory structure should now look like this:  

``` 
# Because directory structures from here on out become much more complicated, 
# I will not be displaying them frequently.
JZ_UCE_RUNS
├── clean-fastq
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
cd clean-fastq/

for i in *;
do
    phyluce_assembly_get_fastq_lengths --input $i/split-adapter-quality-trimmed/ --csv;
done
```
```

```
Each sample averaged about  million reads after illumiprocessing. Not bad. The results are stored in fastqreads.csv.

###<a name="assembly">STEP 2: ASSEMBLE READS AND MATCH PROBES</a>
<a href="https://github.com/faircloth-lab/phyluce">phyluce</a> contains packages for three assembly methods: <a href="http://www.ebi.ac.uk/~zerbino/velvet/">velvet</a>, <a href="http://www.bcgsc.ca/platform/bioinfo/software/abyss">abyss</a>, and <a href="http://trinityrnaseq.sourceforge.net/">trinity</a>. The configuration files you feed into each assembly are similar, so you can make easy modifications to adapt to a different program. I use <a href="http://trinityrnaseq.sourceforge.net/">trinity</a> for this analysis.

I created a config file that directed trinity to wherever my split-adapter-quality-trimmed reads are.

```
# make sure directory paths are correct. otherwise, trinity will not run.
aristochromis-christiae-ID25:/home/jimmyzheng/JZ_UCE_RUNS/clean-fastq/aristochromis-christiae-ID25/split-adapter-quality-trimmed/
astatotilapia-burtoni-ID32:/home/jimmyzheng/JZ_UCE_RUNS/clean-fastq/astatotilapia-burtoni-ID32/split-adapter-quality-trimmed/
aulonocara-stuartgranti-ID26:/home/jimmyzheng/JZ_UCE_RUNS/clean-fastq/aulonocara-stuartgranti-ID26/split-adapter-quality-trimmed/
bathybates-minor-ID90:/home/jimmyzheng/JZ_UCE_RUNS/clean-fastq/bathybates-minor-ID90/split-adapter-quality-trimmed/
boulengerochromis-microlepis-ID79:/home/jimmyzheng/JZ_UCE_RUNS/clean-fastq/boulengerochromis-microlepis-ID79/split-adapter-quality-trimmed/
capidochromis-eucinostomus-ID12:/home/jimmyzheng/JZ_UCE_RUNS/clean-fastq/capidochromis-eucinostomus-ID12/split-adapter-quality-trimmed/
chalinochromis-brichardi-ID57:/home/jimmyzheng/JZ_UCE_RUNS/clean-fastq/chalinochromis-brichardi-ID57/split-adapter-quality-trimmed/
... (continued)
```
I then saved this file in my `JZ_UCE_RUNS` base directory, from which the assembly will be generated. Trinity will create a new directory, which I specified below in the following script. The program runs at a program-defined k-mer value of 25. Smaller k-mer values allows for locating more overlapping sequences but smaller sized subsequences. Thus this can lead to more sequence ambiguities. 25 is considered a small-medium k-mer value, so I also ran trinity with a `min-kmer-coverage` of 60. Higher k-mer values leads to fewer ambiguities, longer contigs, but can take more computing time and space. [INSERT RESULTS HERE]

```
phyluce_assembly_assemblo_trinity \
    --config assembly.conf \
    --output /home/jimmyzheng/JZ_UCE_RUNS/trinity-assemblies \
##  --min-kmer-coverage 60 \ NOTE: Default is 25
    --clean \
    --cores 12
```
Once this finished (took about 4 days), I ran an additional quality check for the assemblies. The results are stored in [contigs.csv](contigs.csv).

```
for i in trinity-assemblies/contigs/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv;
done
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
    --output uce-search-results
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

In my case, <a href="https://www.broadinstitute.org/ftp/pub/assemblies/fish/P_nyererei/PunNye1.0/">*pundamilia nyererei*</a>, as well as the outgroups <a href="ftp://ftp.ensembl.org/pub/release-82/fasta/oreochromis_niloticus/dna/">*oreochromis niloticus*</a> and <a href="https://www.broadinstitute.org/ftp/pub/assemblies/fish/N_brichardi/NeoBri1.0/">*neolamprologus brichardi*</a>, were publicly available assembled genomes, so I just needed to parse and match. Upon web downloading the data, I wrote the following to create an `outgroup-taxon-sets` directory that contains the uce-search-results for just the external data.

```
# to convert the fastas to 2bit files, processable by the next program
faToTwoBit P_nyererei_v1.assembly.fasta p_nyererei.assembly.2bit
faToTwoBit N_brichardi_v1.assembly.fasta n_brichardi.assembly.2bit
faToTwoBit Oreochromis_niloticus.Orenil1.0.dna.toplevel.fa oreochromis-niloticus-genome.2bit

# separate directories for each external dataset
mkdir oreochromis_niloticus_genome; mkdir p_nyererei; mkdir n_brichardi;

mkdir -p outgroup-taxon-sets/outgroup-lastz;

# o_niloticus consists of chromosomal sequences
# p_nyererei and n_brichardi consists of scaffold sequences
phyluce_probe_run_multiple_lastzs_sqlite \
	--db uce-1k-probes.sqlite \
	--output outgroup-lastz \  # output files in subdirectory
    --probefile uce-1k-probes.fasta \
    --chromolist oreochromis_niloticus_genome \
    --scaffoldlist p_nyererei n_brichardi \ # separate by single space
    --genome-base-path ~/JZ_UCE_RUNS/outgroup-taxon-sets \
    --cores 12
```

This will generate DNA alignments for each genome. Then, I sliced the fasta sequences from these lastz files for probe matching. I created a config file for the program to recognize where the original genomes are located.

Config File: [genome-sequence-location.conf] (genome-sequence-location.conf) in base directory

```
[chromos]
oreochromis_niloticus_genome:/home/jimmyzheng/JZ_UCE_RUNS/outgroup-taxon-sets/oreochromis_niloticus_genome/oreochromis_niloticus_genome.2bit

[scaffolds]
p_nyererei:/home/jimmyzheng/JZ_UCE_RUNS/outgroup-taxon-sets/p_nyererei/p_nyererei.2bit
n_brichardi:/home/jimmyzheng/JZ_UCE_RUNS/outgroup-taxon-sets/n_brichardi/n_brichardi.2bit
```

```
cd outgroup-taxon-sets

phyluce_probe_slice_sequence_from_genomes \
	--conf ../genome-sequence-location.conf \
	--lastz outgroup-lastz \
	--output outgroup-fasta \
	--flank=2000 \
	--name-pattern "uce-1k-probes.fasta_v_{}.lastz.clean"
```
Once the fastas for each genome are created, I then matched contigs to probes, which generated the `uce-search-results` subdirectory. Finally, I appended these results to our existing cichlid dataset.

```
cp ../uce-1k-probes.fasta ./
phyluce_assembly_match_contigs_to_probes \
    --contigs outgroup-fasta \
    --probes uce-1k-probes.fasta \
    --output uce-search-results   
```

The data merge requires a config file called [taxon-set.conf] (taxon-set.conf). This is where I indicated which taxa should go into the final matrix generation and pylogenetic analysis.

Config file: [taxon-set.conf] (taxon-set.conf) in base directory

```
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
n_brichardi*
p_nyererei*
oreochromis_niloticus_genome*
# starred species represent external data

```
Script for extending existing cichlid dataset:

```
phyluce_assembly_get_match_counts \
	--locus-db uce-search-results/probe.matches.sqlite \
    --taxon-list-config taxon-set.conf \
    --taxon-group 'extended' \
	--extend outgroup-taxon-sets/uce-search-results/probe.matches.sqlite \
    --output taxon-sets/extended/all-taxa-incomplete.conf \
	--incomplete-matrix 
```

The last step before alignments and trimming is to create an incomplete matrix for **all** of our data, including external genomes. I encountered some issues with naming, as dashes can choke the following program. The safest naming scheme, in my experience, is: `genus_species_###`, with underscores and numbers.

```
phyluce_assembly_get_fastas_from_match_counts \
    --contigs ../../trinity-assemblies/contigs \
    --locus-db ../../uce-search-results/probe.matches.sqlite \
    --extend-locus-contigs ../../outgroup-taxon-sets/outgroup-fasta \
    --extend-locus-db ../../uce-1k-probe-set-extended-outgroups/uce-search-results/probe.matches.sqlite \
    --match-count-output all-taxa-incomplete.conf \
    --output all-taxa-incomplete.fasta \
    --incomplete-matrix all-taxa-incomplete.incomplete \
    --log-path log
```

###<a name="alignment">STEP 4: ALIGNING SEQUENCES AND TRIMMING</a>

Once the incomplete data matrix is created, we need to explode the matrix into fasta files for each species based on each UCE loci sequence. To do this, I ran the following script:

```
cd taxon-sets/extended

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
aristochromis-christiae-ID25.unaligned.fasta,924,644945,697.992424242,4.3301018286,272,1171,706.0,10
astatotilapia-burtoni-ID32.unaligned.fasta,906,573653,633.171081678,4.09758343352,224,1021,642.0,1
aulonocara-stuartgranti-ID26.unaligned.fasta,905,586699,648.286187845,4.26128644469,243,1052,654.0,2
capidochromis-eucinostomus-ID12.unaligned.fasta,899,688248,765.570634038,4.91996539852,266,1226,777.0,45
cheilotilapia-euchilus-ID23.unaligned.fasta,914,633500,693.107221007,4.55016090336,236,1043,702.0,6
... (continued)
```
Now for **aligning the data**. Be sure to specify the right number of taxa, *including external and outgroup taxa*. Otherwise, you might end up creating an extremely incomplete matrix, affecting downstream analyses. Send all log files into the `log` subdirectory.

```
phyluce_align_seqcap_align \
    --fasta all-taxa-incomplete.fasta \
    --output mafft-nexus-edge-trimmed \
    --taxa 37 \
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
2015-11-06 08:13:19,919 - phyluce_align_get_align_summary_data - INFO - ========= Starting phyluce_align_get_align_summary_data =========
2015-11-06 08:13:19,920 - phyluce_align_get_align_summary_data - INFO - Version: 1.5.0
2015-11-06 08:13:19,920 - phyluce_align_get_align_summary_data - INFO - Argument --alignments: /home/jimmyzheng/JZ_UCE_RUNS/taxon-sets/extended/mafft-nexus-internal-trimmed-gblocks
2015-11-06 08:13:19,920 - phyluce_align_get_align_summary_data - INFO - Argument --cores: 12
2015-11-06 08:13:19,920 - phyluce_align_get_align_summary_data - INFO - Argument --input_format: nexus
2015-11-06 08:13:19,920 - phyluce_align_get_align_summary_data - INFO - Argument --log_path: /home/jimmyzheng/JZ_UCE_RUNS/taxon-sets/extended/log
2015-11-06 08:13:19,921 - phyluce_align_get_align_summary_data - INFO - Argument --show_taxon_counts: False
2015-11-06 08:13:19,921 - phyluce_align_get_align_summary_data - INFO - Argument --verbosity: INFO
2015-11-06 08:13:19,921 - phyluce_align_get_align_summary_data - INFO - Getting alignment files
2015-11-06 08:13:19,929 - phyluce_align_get_align_summary_data - INFO - Computing summary statistics using 12 cores
2015-11-06 08:13:22,656 - phyluce_align_get_align_summary_data - INFO - ----------------------- Alignment summary -----------------------
2015-11-06 08:13:22,657 - phyluce_align_get_align_summary_data - INFO - [Alignments] loci:	1,088
2015-11-06 08:13:22,657 - phyluce_align_get_align_summary_data - INFO - [Alignments] length:	650,068
2015-11-06 08:13:22,657 - phyluce_align_get_align_summary_data - INFO - [Alignments] mean:	597.49
2015-11-06 08:13:22,657 - phyluce_align_get_align_summary_data - INFO - [Alignments] 95% CI:	15.10
2015-11-06 08:13:22,658 - phyluce_align_get_align_summary_data - INFO - [Alignments] min:	211
2015-11-06 08:13:22,658 - phyluce_align_get_align_summary_data - INFO - [Alignments] max:	4,101
2015-11-06 08:13:22,661 - phyluce_align_get_align_summary_data - INFO - ------------------------- Taxon summary -------------------------
2015-11-06 08:13:22,661 - phyluce_align_get_align_summary_data - INFO - [Taxa] mean:		31.02
2015-11-06 08:13:22,661 - phyluce_align_get_align_summary_data - INFO - [Taxa] 95% CI:	0.33
2015-11-06 08:13:22,662 - phyluce_align_get_align_summary_data - INFO - [Taxa] min:		3
2015-11-06 08:13:22,662 - phyluce_align_get_align_summary_data - INFO - [Taxa] max:		37
2015-11-06 08:13:22,664 - phyluce_align_get_align_summary_data - INFO - ----------------- Missing data from trim summary ----------------
2015-11-06 08:13:22,664 - phyluce_align_get_align_summary_data - INFO - [Missing] mean:	0.00
2015-11-06 08:13:22,665 - phyluce_align_get_align_summary_data - INFO - [Missing] 95% CI:	0.00
2015-11-06 08:13:22,665 - phyluce_align_get_align_summary_data - INFO - [Missing] min:	0.00
2015-11-06 08:13:22,665 - phyluce_align_get_align_summary_data - INFO - [Missing] max:	0.00
2015-11-06 08:13:22,679 - phyluce_align_get_align_summary_data - INFO - -------------------- Character count summary --------------------
2015-11-06 08:13:22,680 - phyluce_align_get_align_summary_data - INFO - [All characters]	19,721,992
2015-11-06 08:13:22,680 - phyluce_align_get_align_summary_data - INFO - [Nucleotides]		19,415,735
2015-11-06 08:13:22,682 - phyluce_align_get_align_summary_data - INFO - ---------------- Data matrix completeness summary ---------------
2015-11-06 08:13:22,682 - phyluce_align_get_align_summary_data - INFO - [Matrix 50%]		1054 alignments
2015-11-06 08:13:22,682 - phyluce_align_get_align_summary_data - INFO - [Matrix 55%]		1050 alignments
2015-11-06 08:13:22,682 - phyluce_align_get_align_summary_data - INFO - [Matrix 60%]		1036 alignments
2015-11-06 08:13:22,682 - phyluce_align_get_align_summary_data - INFO - [Matrix 65%]		1021 alignments
2015-11-06 08:13:22,682 - phyluce_align_get_align_summary_data - INFO - [Matrix 70%]		984 alignments
2015-11-06 08:13:22,683 - phyluce_align_get_align_summary_data - INFO - [Matrix 75%]		930 alignments
2015-11-06 08:13:22,683 - phyluce_align_get_align_summary_data - INFO - [Matrix 80%]		835 alignments
2015-11-06 08:13:22,683 - phyluce_align_get_align_summary_data - INFO - [Matrix 85%]		649 alignments
2015-11-06 08:13:22,683 - phyluce_align_get_align_summary_data - INFO - [Matrix 90%]		514 alignments
2015-11-06 08:13:22,683 - phyluce_align_get_align_summary_data - INFO - [Matrix 95%]		240 alignments
2015-11-06 08:13:22,683 - phyluce_align_get_align_summary_data - INFO - ------------------------ Character counts -----------------------
2015-11-06 08:13:22,683 - phyluce_align_get_align_summary_data - INFO - [Characters] '-' is present 306,257 times
2015-11-06 08:13:22,684 - phyluce_align_get_align_summary_data - INFO - [Characters] 'A' is present 5,424,326 times
2015-11-06 08:13:22,684 - phyluce_align_get_align_summary_data - INFO - [Characters] 'C' is present 4,253,925 times
2015-11-06 08:13:22,684 - phyluce_align_get_align_summary_data - INFO - [Characters] 'G' is present 4,303,860 times
2015-11-06 08:13:22,684 - phyluce_align_get_align_summary_data - INFO - [Characters] 'T' is present 5,433,624 times
2015-11-06 08:13:22,684 - phyluce_align_get_align_summary_data - INFO - ========= Completed phyluce_align_get_align_summary_data ========
```
Our 75% matrix looks good at 930 alignments. We can also choose not to trim and see how our alignments look. In our case, they do not look much different, so we move forward. The script to run with just internal trimming is as follows:

```
phyluce_align_seqcap_align \
    --fasta all-taxa-incomplete.fasta \
    --output mafft-nexus-internal-trimmed \
    --taxa 37 \
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
    --taxa 37 \ # make sure this number has been accurate throughout
    --percent 0.75 \ # to specify 75% incompleteness
    --output mafft-nexus-internal-trimmed-gblocks-clean-75p \
    --cores 12 \
    --log-path log

# FOR RAXML    
phyluce_align_format_nexus_files_for_raxml \
    --alignments mafft-nexus-internal-trimmed-gblocks-clean-75p \
    --output mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml \
    --charsets \
    --log-path log

# SAME ALIGNMENT FOR EXABAYES
cp -R mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml mafft-nexus-internal-trimmed-gblocks-clean-75p-exabayes
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
│	     ├── mafft-nexus-internal-trimmed-gblocks-clean-75p  
│	     ├── mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml  
│	     ├── mafft-nexus-internal-trimmed-gblocks-clean-75p-exabayes  
│	     └── log
└── ... (continued)
```

Then, I began my RAxML analyses. 

```
raxmlHPC-PTHREADS-SSE3 \
    -m GTRGAMMA \
    -N 100 \
    -p 25365 \
    -n best \
    -s mafft-nexus-internal-trimmed-gblocks-clean-75p.phylip \
    -T 12 \
    -o simochromis_babaulti_ID48
    
# SOME PROBLEMS:
# -too many reps
# -did not specify partitions (explained below)
# -outgroup not distant enough
```
At first, I made the mistake of not including the partitions. Because there are so many UCE loci, it is incredibly important to specify where each loci begins and ends, because without those data, RAxML treats the whole sequence for each taxa as a single loci. This confuses phylogenetic inference. The resultant tree can be found at [RAxML_extended.cichlid.tre] (RAxML_extended.cichlid.tre). As you can see in FigTree, the bootstrap support values sucked.

To get the right partitions, I went into the file `mafft-nexus-internal-trimmed-gblocks-clean-75p.charsets` in the `mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml` directory. It gave all the partition information that I needed:

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
Copying this into another file called `aln.part`, I sorted by the UCE loci positions in increasing order and modified the format.

```
DNA, uce-655=1-587
DNA, uce-657=588-1108
DNA, uce-1027=1109-1652
DNA, uce-279=1653-2321
DNA, uce-620=2322-2858
DNA, uce-854=2859-3563
DNA, uce-986=3564-4035
DNA, uce-661=4036-4610
DNA, uce-1115=4611-5157
DNA, uce-73=5158-5759
... (continued)
```

Then I reran the RAxML analysis in a separate folder called `mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml-mod`.

I am currently in the process of running the new tree...

# TO BE CONTINUED
