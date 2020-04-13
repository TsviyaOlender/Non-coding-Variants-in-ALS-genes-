# Non-coding-Variants-in-ALS-genes-
Scripts and data files used for the annotation of 3'UTR and microRNA variants

#
The perl script pipeline_joint_V2.pl performs a detailed annotation of VCF file.
usage:
perl pipeline_joint_V2.pl VCFfile cohortName
The script is based havily on annovar (https://doc-openbio.readthedocs.io/projects/annovar/en/latest/). Annovar scripts have to be placed under folder "annovar".
The following files have to be downloaded from annovar site and to be placed under annovar/humandb:
refGene,avsnp144,1000g2014oct_eur,esp6500siv2_ea,cg69,phastConsElements46way,genomicSuperDups,dbnsfp30a,wgEncodeRegDnaseClustered,
wgEncodeRegTfbsClustered,clinvar_20150629,exac03,

Our special annotation for microRNA and 3UTR variants: 
targetScanS,preMiRNA,matureMiRNA,seedMiRNA,TargetS16_conSites_poorConsFam,TargetS16_nonconsvSites_consvFam
TargetS16_noncosvSites_broadConsFam,TargetS16_sites,GenCode3UTR
The files PITA_Ts_run and PITA_top_f- include prediction of miRNA binding site loss using PITA, performed on a restricted list of ALS related genes.
GenCode3UTR- annotates specifically the 3UTRs of ALS related genes.
The script also scans a file called eXome_db which contains variations of the Israeli population. The online version contains dummy information due to ethical reasons.

#
Additional scripts:
GenotypeToSummary.pl
Finding which individuals carry specific mutations

filter_annotation_file.pl
Removes multi-allelic variants from the annotation file
