#!/usr/bin/perl/ -w
$|++;
use strict;
#use Math::Random; #---For generating normal distribution random numbers
use File::Path; #---for removing tmp GNUPLOT dat files
use Time::HiRes qw (sleep);
use List::Util qw ( sum );

######################################################################################################################################################
#
#	Description
#		This script is mainly built based on geneModelBuilder. This is a perl script to reconstruct the alternatively splicing isoforms of bona fide reference gene models. 
#	It will not reconstruct new transcripts outside the reference gene models but instead look for alternative splicing junction WITHIN the bona fide gene models.
#
#	splicing;
#
#	Input
#		--refGffPath=			path of the reference GFF, usually generated from geneModelReviseInfoCombinner, the bona fide gff;
#		--canJunctionBedPath=	the bed of the canoical junction; normally generated from HMMSplicerBEDToSAMParser.pl;
#		--nonCanJunctionBedPath=the bed of the non-canonial cannoical junction; normally generated from HMMSplicerBEDToSAMParser.pl;
#		--pileupPath=			the pileup file generated from pileup counter, used to determine the ration of splice and unspliced read flanking the splicing jucntions;
#		--minSplicingRatio=		the minimum ratio of spliced to unspliced read for a junction to be regarded as a "major" junction;
#		--minMajorSupportRead=	the minimum number of supporting read for a junction to be regarded as a "major" junction;
#		--intronBoundCovPath=	"no" or auto or the path of a file that contains the intron bound coverage of the junctions. If "no", the coverage will be calculated from the pileu file; If "auto", the script will search for it presence in the default path;
#		--geneticCodePath=		the path to the genetic code table;
#		--cntgFastaPath=		the path of the cntg Fasta;
#		--geneticCodeID=		the ID of the genetic code model; default = 1, i.e standard code
#		--ATGOnly=				"yes" or "no", if "yes", only ATG will be treated as start codon and the alternative start codon specified in the codon table will not be considered; default = yes;
#		--unqPileupPath=		the pileup file generated from pile up counter for gSimulation dateset, use for calculating the uniqueness; use "no" to use precalculated ranges;
#		--precalUnqLoCmPath=		the path of the precalculated uniqueness; will be ignored if unqPileupPath is not equal to "no";
#		--boundWidth			the width of the bound that is used to calculate uniqueness and get flanking sequence around the junctions; default = 50;
#		--missModJunctComPath=	path of a file that contains the missed and modified junction comments;
#		--ctrgyToCount=			ctgry in Gff file for counting; default = all;
#		--countFturCovInfoPath=	path of Cov per Nt info of all countFutr, for calculating the splicing efficiency; use "no" to read from pileupPath
#		--pctCutLim=			the low percentile of all refHit splicing efficiency to be used as the cutoff to define "often" or "rare"; default = 0.9
#		--strndSpecificPileup	"yes" or "no"; if "yes", only the strand specific cov will be taken into account; default = no
#		--lowComRegionPath		the path of the output of dustmakser; use "no" to skip;
#		--outDir= 				output directory
#
#	Output
#		
#
#	Usage
#		perl alternativeSplicer_v0.2.pl --refGffPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pipeLines/marathon/resources/genome/EHI_v13_original.gff --junctionBedPath=/Volumes/A_MPro2TB/NGS/full/pooledData/HMMSplicer/HMMSplicerBEDToSAMParser/HMMSplicer/junctionInfo/HMMSplicer.filter.can.unique.dup.bed --outDir=./EHI_v13Ref/ --cntgFastaPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pipeLines/marathon/resources/genome/EHI_v13.fa --geneticCodePath=/Volumes/CCHON1TBA/otherResources/NCBIGeneticCode.txt --ATGOnly=yes --geneticCodeID=1 --boundWidth=50 --missModJunctComPath=missModJunctComments.txt --minSplicingRatio=0.3 --minMajorSupportRead=3 --GffForCountingPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/geneModelReviseInfoCombinner/test/all.revised.gff --pileupPath=no --unqPileupPath=no --intronBoundCovPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/alternativeSplicer/EHI_v13Ref/preCalInfo/intronBoundCov.txt --countFturCovInfoPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/alternativeSplicer/EHI_v13Ref/preCalInfo/countFturCovPerNt.txt --precalUnqLoCmPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/alternativeSplicer/EHI_v13Ref/preCalInfo/precalculated.unq.txt --pctCutLim=0.95
#		perl alternativeSplicer_v0.2.pl --refGffPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/geneModelReviseInfoCombinner/test/all.revised.mRNA.gff --junctionBedPath=/Volumes/A_MPro2TB/NGS/full/pooledData/HMMSplicer/HMMSplicerBEDToSAMParser/HMMSplicer/junctionInfo/HMMSplicer.filter.can.unique.dup.bed --outDir=./allRevisedmRNA/ --cntgFastaPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pipeLines/marathon/resources/genome/EHI_v13.fa --geneticCodePath=/Volumes/CCHON1TBA/otherResources/NCBIGeneticCode.txt --ATGOnly=yes --geneticCodeID=1 --boundWidth=50 --missModJunctComPath=missModJunctComments.txt --minSplicingRatio=0.3 --minMajorSupportRead=3 --GffForCountingPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/geneModelReviseInfoCombinner/test/all.revised.gff --pileupPath=no --unqPileupPath=no --intronBoundCovPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/alternativeSplicer/allRevisedmRNA/preCalInfo/intronBoundCov.txt --countFturCovInfoPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/alternativeSplicer/allRevisedmRNA/preCalInfo/countFturCovPerNt.txt --precalUnqLoCmPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/alternativeSplicer/allRevisedmRNA/preCalInfo/precalculated.unq.txt --pctCutLim=0.95
#		perl alternativeSplicer_v0.2.pl --refGffPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/geneModelReviseInfoCombinner/test/unq.bonaFide.revised.gff --junctionBedPath=/Volumes/A_MPro2TB/NGS/full/pooledData/HMMSplicer/HMMSplicerBEDToSAMParser/HMMSplicer/junctionInfo/HMMSplicer.filter.can.unique.dup.bed --outDir=./bfmRNAAltSplice/ --cntgFastaPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pipeLines/marathon/resources/genome/EHI_v13.fa --geneticCodePath=/Volumes/CCHON1TBA/otherResources/NCBIGeneticCode.txt --ATGOnly=yes --geneticCodeID=1 --boundWidth=50 --missModJunctComPath=missModJunctComments.txt --minSplicingRatio=0.3 --minMajorSupportRead=3 --GffForCountingPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/geneModelReviseInfoCombinner/test/all.revised.gff --pileupPath=no --unqPileupPath=no --intronBoundCovPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/alternativeSplicer/bfmRNAAltSplice/preCalInfo/intronBoundCov.txt --countFturCovInfoPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/alternativeSplicer/bfmRNAAltSplice/preCalInfo/countFturCovPerNt.txt --precalUnqLoCmPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/transcriptModel/alternativeSplicer/bfmRNAAltSplice/preCalInfo/precalculated.unq.txt --pctCutLim=0.95
#
#	Version
#		v0.2
#		-Added function to compare the ref and the NGS juncts
#
#		v0.3
#		--corrected the method to construct alternative transcripts, instead of sort all boundaries.
#		--added strndSpecificPileup option
#
#		v0.4
#		--added analyses for non-cannonical junction
#		--added quantitave XY weblog entropy vs splicing efficiency
#		--will use dreme to look for motifs in surrounding exon
#
#		v0.5
#		--will output geneBasedAntisense.info.txt;
#		--added function to detect splicing efficieny vs proportion of preserving ORF;
#		--added mean number of junctions within a cluster versus representing junction supporting read number interval
#		--will output the abundance of exon skipped and non exon skipped multiple intron genes
#
#		v0.6
#		--witin sub calculateMeanJunctionNumWithinClusterRepJunctRdNumInterval., added estimation of splicing error rate and try to pick up the non-stochastic alternative splicing junctions;
#		--try to pickup the clusters that are likely to be purely noise based on consensuse value of the alternative junctions.
#		--will calculate whether the junction is 3n, 3n+1, 3n+2, reported in all_NGS_junct_log.txt.
#		--added function to investigate whether high efficiency altSites are more likely to occurr on multiple introns of the same transcript
#		--will search for seq motifs in non-Stochastic altSite and intron create in both dreme and weblogo, in dremeNonStochasticJunctions function;
#
######################################################################################################################################################

#==========================================================Main body starts==========================================================================#

#----------Read the parameters----------#
use vars qw ($dummy $refGffPath $canJunctionBedPath $nonCanJunctionBedPath $pileupPath $minSplicingRatio $intronBoundCovPath $minMajorSupportRead $geneticCodePath $cntgFastaPath $geneticCodeID $ATGOnly $unqPileupPath $boundWidth $missModJunctComPath $precalUnqLoCmPath $countFturCovInfoPath $pctCutLim $strndSpecificPileup $lowComRegionPath $ctrgyToCountHsh_ref $outDir);
my ($refGffPath, $canJunctionBedPath, $nonCanJunctionBedPath, $pileupPath, $minSplicingRatio, $intronBoundCovPath, $minMajorSupportRead, $geneticCodePath, $cntgFastaPath, $geneticCodeID, $ATGOnly, $unqPileupPath, $boundWidth, $missModJunctComPath, $precalUnqLoCmPath, $countFturCovInfoPath, $pctCutLim, $strndSpecificPileup, $lowComRegionPath, $ctrgyToCountHsh_ref, $outDir) = readParameters();
printCMDLogOrFinishMessage("CMDLog");

#----read the seq
my $fastaHsh_ref = readMultiFasta($cntgFastaPath);

#----get the uniquess
my $unqLoComByCntgPosHsh_ref = getUnqAndLoComInfo($fastaHsh_ref, $unqPileupPath, $lowComRegionPath);

#----------Read the ref Gff
my ($refRngSSHsh_ref, $refCtgryByGeneHsh_ref, $refStrndHsh_ref, $refRngXSHsh_ref, $refExonRngByGeneHsh_ref, $refIntronRngByGeneHsh_ref, $refJunctStrByIntronIDHsh_ref, $refIntronRngXS_ref, $refStrndByIntronHsh_ref, $refNameByGeneHsh_ref, $refGeneIDByJunctStrHsh_ref, $refGeneCntgByGeneHsh_ref) = readGff($refGffPath, $ctrgyToCountHsh_ref);

#----------Read the Gff for category counting
my ($fturCountRngSSHsh_ref, $fturCountCtgryByGeneHsh_ref, $fturCountStrndHsh_ref, $fturCountRngXSHsh_ref, $fturCountExonRngByGeneHsh_ref, $fturCountIntronRngByGeneHsh_ref, $fturCountJunctStrByIntronIDHsh_ref, $fturCountIntronRngXS_ref, $fturCountStrndByIntronHsh_ref, $fturCountNameByGeneHsh_ref, $fturCountGeneIDByJunctStrHsh_ref, $fturCountGeneCntgByGeneHsh_ref) = readGff($refGffPath, $ctrgyToCountHsh_ref);

my $intronlessGeneCodingSeqHsh_ref = getCodingSequence($fastaHsh_ref, $fturCountStrndHsh_ref, $fturCountExonRngByGeneHsh_ref, $fturCountGeneCntgByGeneHsh_ref, $fturCountCtgryByGeneHsh_ref);

#----------Read the BED
my ($NGSJunctStrndHsh_ref, $junctRngHsh_ref, $NGSJunctCntgHsh_ref, $junctScoreHsh_ref, $junctReadNumHsh_ref, $junctIntronRngHsh_ref, $NGSJunctBEDCoverRngHsh_ref) = readJunctBED($canJunctionBedPath);

#----------Read the non-cannonical BED
my ($nonCanNGSJunctStrndHsh_ref, $nonCanJunctRngHsh_ref, $nonCanNGSJunctCntgHsh_ref, $nonCanJunctScoreHsh_ref, $nonCanJunctReadNumHsh_ref, $nonCanJunctIntronRngHsh_ref, $nonCanNGSJunctBEDCoverRngHsh_ref) = readJunctBED($nonCanJunctionBedPath);

#---------Compute the overlapping canonical and non-canonical junctions
my ($dummy11, $dummy12, $XSOvrlpNonCanNGSJCanNGSJHsh_ref, $XSOvrlpCanNGSJNonCanNGSJHsh_ref) = findFturOverlap($nonCanJunctIntronRngHsh_ref, $junctIntronRngHsh_ref, $nonCanNGSJunctStrndHsh_ref, $NGSJunctStrndHsh_ref,  "y", "y", "all", "all", "Scanning the overlapping between cannoical and noncannoical NGS junctions");

#---------summarize the overlapping between canonical and noncanoical junctions
my ($nonCanNGSJunctInfoHsh_ref, $nonCanNGSJunctFilterStrndHsh_ref) = summarizeNonCanAndCanNGSSOverlapping($XSOvrlpNonCanNGSJCanNGSJHsh_ref, $nonCanNGSJunctStrndHsh_ref, $NGSJunctStrndHsh_ref, $nonCanNGSJunctCntgHsh_ref, $nonCanJunctScoreHsh_ref, $nonCanJunctReadNumHsh_ref);
my %nonCanNGSJunctInfoHsh = %{$nonCanNGSJunctInfoHsh_ref};

#---------Compute the overlapping junctions to define overlapping whole BED ranges (implemented in v0.2)
my ($SSOvrlpJunctBEDRngHsh_ref, $dummy1, $XSOvrlpJunctBEDRngHsh_ref, $dummy2) = findFturOverlap($junctRngHsh_ref, $junctRngHsh_ref, $NGSJunctStrndHsh_ref, $NGSJunctStrndHsh_ref, "n", "y", "all", "all", "Scanning self-overlapping of all BED ranges");

#-------define the overlapping junct BED ranges that is used for clustering with the transfrags
my ($junctBEDRngClusterJunctHsh_ref, $junctBEDRngClusterNameByJunctHsh_ref, $junctBEDClusterRngHsh_ref, $junctBEDClusterStrndHsh_ref, $junctBEDClusterCntgHsh_ref) = defineOverlappingJunctBEDRng($SSOvrlpJunctBEDRngHsh_ref, $junctRngHsh_ref, $NGSJunctCntgHsh_ref, $NGSJunctStrndHsh_ref);

#---------Compute the overlapping junctions (only intron but not the whole range)
my ($SSOvrlpJunctIntronHsh_ref, $dummy3, $XSOvrlpJunctHsh_ref, $dummy4) = findFturOverlap($junctIntronRngHsh_ref, $junctIntronRngHsh_ref, $NGSJunctStrndHsh_ref, $NGSJunctStrndHsh_ref, "n", "y", "all", "all", "Scanning the self-overlapping of introns");

#---------Compute the overlapping junctions with the reference junctions
my ($SSOvrlpRefJNGSJHsh_ref, $SSOvrlpNGSJRefJHsh_ref, $XSOvrlpRefJNGSJHsh_ref, $XSOvrlpNGSJRefJHsh_ref) = findFturOverlap($refIntronRngXS_ref, $junctIntronRngHsh_ref, $refStrndByIntronHsh_ref, $NGSJunctStrndHsh_ref,  "y", "y", "all", "all", "Scanning the overlapping between NGS junctions with the reference junctions");

#---------Compute the overlapping nonCanonical junctions with the reference junctions
my ($SSOvrlpRefJNonCanNGSJHsh_ref, $SSOvrlpNonCanNGSJRefJHsh_ref, $XSOvrlpRefJNonCanNGSJHsh_ref, $XSOvrlpNonCanNGSJRefJHsh_ref) = findFturOverlap($refIntronRngXS_ref, $nonCanJunctIntronRngHsh_ref, $refStrndByIntronHsh_ref, $nonCanNGSJunctStrndHsh_ref,  "y", "y", "all", "all", "Scanning the overlapping between nonCanonical NGS junctions with the reference junctions");

#---------Compute the overlapping junctions with the counting features
my ($SSOvrlpNGSJFturCountHsh_ref, $dummy7a, $XSOvrlpNGSJFturCountHsh_ref, $dummy7b) = findFturOverlap($junctIntronRngHsh_ref, $fturCountRngXSHsh_ref, $NGSJunctStrndHsh_ref, $fturCountStrndHsh_ref, "y", "y", "all", "all", "Scanning the overlapping between NGS junctions with the counting features");

#---------Compute the overlapping junctions with the counting features
my ($SSOvrlpNonCanNGSJFturCountHsh_ref, $dummy8a, $XSOvrlpNonCanNGSJFturCountHsh_ref, $dummy8b) = findFturOverlap($nonCanJunctIntronRngHsh_ref, $fturCountRngXSHsh_ref, $nonCanNGSJunctStrndHsh_ref, $fturCountStrndHsh_ref, "y", "y", "all", "all", "Scanning the overlapping between nonCanonical NGS junctions with the counting features");

#----------get Junctionn Unq And Seq of cannonical junctions
my ($NGSJunctUnqHsh_ref, $NGSJunctSeqHsh_ref, $refJunctInfoHsh_ref);
my $exonSeqLength = 10;
($NGSJunctUnqHsh_ref, $NGSJunctSeqHsh_ref, $refJunctInfoHsh_ref, $nonCanNGSJunctInfoHsh_ref) = getJunctUnqAndSeq($NGSJunctCntgHsh_ref, $unqLoComByCntgPosHsh_ref, $NGSJunctStrndHsh_ref, $fastaHsh_ref, $refGeneIDByJunctStrHsh_ref, $refStrndHsh_ref, $SSOvrlpRefJNGSJHsh_ref, \%nonCanNGSJunctInfoHsh, $exonSeqLength);
%nonCanNGSJunctInfoHsh = %{$nonCanNGSJunctInfoHsh_ref};

#---------Compute the overlapping reference and canonical junctions
my ($SSOvrlpRefTrnscptNGSJHsh_ref, $dummy9a, $XSOvrlpRefTrnscptNGSJHsh_ref, $dummy10a) = findFturOverlap($refRngXSHsh_ref, $junctIntronRngHsh_ref, $refStrndHsh_ref, $NGSJunctStrndHsh_ref, "y", "n", "all", "all", "Scanning overlapping between reference and Cannonical NGS junctions ");

#---------Compute the overlapping reference and nonCanonical junctions
my ($SSOvrlpRefTrnscptNonCanNGSJHsh_ref, $SSOvrlpNonCanNGSJRefTrnscptHsh_ref, $XSOvrlpRefTrnscptNonCanNGSJHsh_ref, $XSOvrlpNonCanNGSJRefTrnscptHsh_ref) = findFturOverlap($refRngXSHsh_ref, $nonCanJunctIntronRngHsh_ref, $refStrndHsh_ref, $nonCanNGSJunctStrndHsh_ref, "y", "n", "all", "all", "Scanning overlapping between reference and nonCanonical NGS junctions ");

#---------scanJunctionIntronBoundCoverage
my ($junctIntronBoundCovHsh_ref, $countFturCovInfoHsh_ref) = scanJunctionIntronBoundCoverage($junctIntronRngHsh_ref, $fturCountRngXSHsh_ref, $fturCountExonRngByGeneHsh_ref);

#----------define the intron clusters
my ($junctSplicingRatioHsh_ref, $majorJunctAvgSplcingRatioHsh_ref, $superJunctHsh_ref, $clusterAllJunctHsh_ref, $jClusterNameByJunctHsh_ref, $jClusterInfoHsh_ref, $jClusterSSOverlapHsh_ref, $totalExonSkippingTypeHsh_ref, $majorExonSkippingTypeHsh_ref, $exactExonSkippingClusterHsh_ref, $splicingSiteDiffAllHsh_ref, $splicingSiteDiffHshByClusterHsh_ref, $junctOnProminentIsofmHsh_ref, $superJunctOvrlpClusterHsh_ref, $NGSJunctInfoHsh_ref) = defineIntronClusters($SSOvrlpJunctIntronHsh_ref, $NGSJunctStrndHsh_ref, $junctScoreHsh_ref, $junctReadNumHsh_ref, $junctIntronBoundCovHsh_ref, $junctIntronRngHsh_ref, $XSOvrlpJunctHsh_ref, $refGeneIDByJunctStrHsh_ref, $NGSJunctUnqHsh_ref, $NGSJunctSeqHsh_ref, $NGSJunctBEDCoverRngHsh_ref);

#----count the NGS JunctStr overlap
my ($oftenSplicingEffCutoff, $geneBasedJunctCountHsh_ref);
($NGSJunctInfoHsh_ref, $oftenSplicingEffCutoff, $geneBasedJunctCountHsh_ref) = countNGSJunctOverlapFtur($NGSJunctInfoHsh_ref, $SSOvrlpNGSJFturCountHsh_ref, $XSOvrlpNGSJFturCountHsh_ref, $fturCountCtgryByGeneHsh_ref, $countFturCovInfoHsh_ref, $junctReadNumHsh_ref, $SSOvrlpJunctIntronHsh_ref);
my %NGSJunctInfoHsh = %{$NGSJunctInfoHsh_ref};

printGeneBasedJunctCount($geneBasedJunctCountHsh_ref, $refNameByGeneHsh_ref);

#----count the noncanoical NGS JunctStr overlap
$nonCanNGSJunctInfoHsh_ref = summarizeNonCanNGSJunctOverlapFtur(\%nonCanNGSJunctInfoHsh, $SSOvrlpNonCanNGSJFturCountHsh_ref, $XSOvrlpNonCanNGSJFturCountHsh_ref, $countFturCovInfoHsh_ref, $nonCanJunctReadNumHsh_ref, $oftenSplicingEffCutoff);
%nonCanNGSJunctInfoHsh = %{$nonCanNGSJunctInfoHsh_ref};

#----read the codon table
my ($selectedIDCodonTableHsh_ref, $selectedIDStartCodonHsh_ref) = getCodonTable($geneticCodePath, $geneticCodeID);

#----summarize Alternative Splicing On RefTranscript
my ($altTrnscptInfoHsh_ref, $altJunctStrInfoHsh_ref);
($altTrnscptInfoHsh_ref, $altJunctStrInfoHsh_ref, $refJunctInfoHsh_ref, $NGSJunctInfoHsh_ref) = summarizeAlternativeSplicingOnRefTranscript($refIntronRngXS_ref, $refRngXSHsh_ref, $refStrndHsh_ref, $junctScoreHsh_ref, $junctReadNumHsh_ref, $NGSJunctCntgHsh_ref, $SSOvrlpRefJNGSJHsh_ref, $XSOvrlpRefJNGSJHsh_ref, $SSOvrlpRefTrnscptNGSJHsh_ref, $XSOvrlpRefTrnscptNGSJHsh_ref, $superJunctHsh_ref, $clusterAllJunctHsh_ref, $jClusterNameByJunctHsh_ref, $jClusterInfoHsh_ref, $totalExonSkippingTypeHsh_ref, $superJunctOvrlpClusterHsh_ref, $SSOvrlpJunctIntronHsh_ref, $refExonRngByGeneHsh_ref, $refJunctStrByIntronIDHsh_ref, $refIntronRngByGeneHsh_ref, $jClusterSSOverlapHsh_ref, $selectedIDCodonTableHsh_ref, $selectedIDStartCodonHsh_ref, $fastaHsh_ref, $refNameByGeneHsh_ref, $refCtgryByGeneHsh_ref, $countFturCovInfoHsh_ref, $NGSJunctInfoHsh_ref, $SSOvrlpNGSJRefJHsh_ref, $refJunctInfoHsh_ref);
my %refJunctInfoHsh = %{$refJunctInfoHsh_ref};
%NGSJunctInfoHsh = %{$NGSJunctInfoHsh_ref};

my $maxUnq = 1;
my $minScore = 1050;
my $minReadNum = 1;
my $maxLoCmPrptn = 0.1;
my $filternonCanNGSJunctFilterStrndHsh_ref = filterAndPrintNonCanNGSJInfo(\%nonCanNGSJunctInfoHsh, $maxUnq, $minScore, $minReadNum, $maxLoCmPrptn);

$maxUnq = 1;
$minScore = 1050;
$minReadNum = 1;
$maxLoCmPrptn = 0.1;
my $maxPolyBasePrptn = 0.7;
summarizeBothCanNonCanNGSJunctAltSiteOnRefJunct($SSOvrlpNonCanNGSJRefTrnscptHsh_ref, $SSOvrlpRefJNGSJHsh_ref, $SSOvrlpRefJNonCanNGSJHsh_ref, \%NGSJunctInfoHsh, \%refJunctInfoHsh, \%nonCanNGSJunctInfoHsh, $maxUnq, $maxLoCmPrptn, $maxPolyBasePrptn);

#------readRefJunctComment
my $refJunctCommentHsh_ref = readRefJunctComment();

#-------plotSplicingEfficiencyVersusSiteEntropy
$NGSJunctInfoHsh_ref = plotSplicingEfficiencyVersusSiteEntropy(\%NGSJunctInfoHsh, 50, $intronlessGeneCodingSeqHsh_ref);
%NGSJunctInfoHsh = %{$NGSJunctInfoHsh_ref};

#---splicing efficieny and the proportion of full ORF
my $numJunctInterval = 50;
my $minReadNumForPrprtnORF = 3;
calculateProportionOfCodingAtSplicingEfficienyInterval($NGSJunctInfoHsh_ref, $numJunctInterval, $minReadNumForPrprtnORF);

#---jucntion num within cluster and representative read num
my $intervalBin = 2;
my $logScale = 2;
calculateMeanJunctionNumWithinClusterRepJunctRdNumInterval($jClusterInfoHsh_ref, $intervalBin, $logScale, $NGSJunctInfoHsh_ref);

#---correlation between intron creation and transcript abundance
calculateProportionOfIntronCreationAtAbundanceInterval($NGSJunctInfoHsh_ref, $countFturCovInfoHsh_ref, $refStrndHsh_ref, $intervalBin, $logScale);

#---disabled temporarily, it causes crash
my $minPrmntReadNumRatioAsNonStochasticAltSite = 0.13;
my $minSenseSplicingEfficiencyAsNonStochasticIntronCreate = 0.02;
dremeNonStochasticJunctions(\%NGSJunctInfoHsh, $intronlessGeneCodingSeqHsh_ref, $minPrmntReadNumRatioAsNonStochasticAltSite, $minSenseSplicingEfficiencyAsNonStochasticIntronCreate);

statisticsOfMultiIntronTranscripts($refGeneIDByJunctStrHsh_ref, $altTrnscptInfoHsh_ref, $altJunctStrInfoHsh_ref, $junctReadNumHsh_ref, $countFturCovInfoHsh_ref, $NGSJunctInfoHsh_ref, $jClusterInfoHsh_ref, $minPrmntReadNumRatioAsNonStochasticAltSite);

#-------sortAndPrintAllJunctInfo
my $allSeqPathByJunctTypeHsh_ref = sortAndPrintAllJunctInfo(\%NGSJunctInfoHsh, \%refJunctInfoHsh, $refJunctCommentHsh_ref, $jClusterInfoHsh_ref);

calculateEntropyAndPlotLogo($allSeqPathByJunctTypeHsh_ref);

#------print hit log----------#
printCMDLogOrFinishMessage("finishMessage");

exit;
#========================================================= Main body ends ===========================================================================#

########################################################################## readParameters
sub readParameters {
	
	$outDir = "./";
	$pctCutLim = 0.9;
	$missModJunctComPath = "no";
	$intronBoundCovPath = "auto";
	$countFturCovInfoPath = "auto";
	$precalUnqLoCmPath = "auto";
	$pileupPath = "no";
	$unqPileupPath = "no";
	$strndSpecificPileup = "no";
	$boundWidth = 50;
	
	my %ctrgyToCountHsh;
	
	foreach my $param (@ARGV) {

		if ($param =~ m/--refGffPath=/) {$refGffPath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--canJunctionBedPath=/) {$canJunctionBedPath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--nonCanJunctionBedPath=/) {$nonCanJunctionBedPath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--pileupPath=/) {$pileupPath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--minSplicingRatio=/) {$minSplicingRatio = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--intronBoundCovPath=/) {$intronBoundCovPath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--minMajorSupportRead=/) {$minMajorSupportRead = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--geneticCodePath=/) {$geneticCodePath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--cntgFastaPath=/) {$cntgFastaPath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--geneticCodeID=/) {$geneticCodeID = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--ATGOnly=/) {$ATGOnly = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--unqPileupPath=/) {$unqPileupPath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--boundWidth=/) {$boundWidth = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--missModJunctComPath=/) {$missModJunctComPath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--ctrgyToCount=/) {my $ctrgyToCount = substr ($param, index ($param, "=")+1); $ctrgyToCountHsh{$ctrgyToCount}++;}
		elsif ($param =~ m/--precalUnqLoCmPath=/) {$precalUnqLoCmPath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--countFturCovInfoPath=/) {$countFturCovInfoPath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--pctCutLim=/) {$pctCutLim = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--strndSpecificPileup=/) {$strndSpecificPileup = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--lowComRegionPath=/) {$lowComRegionPath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--outDir=/) {$outDir = substr ($param, index ($param, "=")+1);}
	}
	
	$ctrgyToCountHsh{"all"}++ if ((keys %ctrgyToCountHsh) == 0);
	
	#---have to read both other may may cause problem
	if ((($lowComRegionPath ne "no") and ($unqPileupPath eq "no")) or (($lowComRegionPath eq "no") and ($unqPileupPath ne "no"))) {
		die "Both lowComRegionPath and unqPileupPath have to be no or specify a path\n";
	}
	
	open (TEST, "$refGffPath") || die "Cannot open $refGffPath: $!\n"; close TEST;
	open (TEST, "$canJunctionBedPath") || die "Cannot open $canJunctionBedPath: $!\n"; close TEST;
	open (TEST, "$nonCanJunctionBedPath") || die "Cannot open $nonCanJunctionBedPath: $!\n"; close TEST;
	if ($lowComRegionPath ne "no") {
		open (TEST, "$lowComRegionPath") || die "Cannot open $lowComRegionPath: $!\n"; close TEST;
	}
	
	open (TEST, "$geneticCodePath") || die "Cannot open $geneticCodePath: $!\n"; close TEST;
	open (TEST, "$cntgFastaPath") || die "Cannot open $cntgFastaPath: $!\n"; close TEST;
	if ($missModJunctComPath ne "no") {
		open (TEST, "$missModJunctComPath") || die "Cannot open $missModJunctComPath: $!\n"; close TEST;
	}

	if ($pileupPath ne "no") {
		open (TEST, "$pileupPath") || die "Cannot open $pileupPath: $!\n"; close TEST;	
	} else {
		$intronBoundCovPath = "$outDir/preCalInfo/intronBoundCov.txt" if ($intronBoundCovPath eq "auto");
		$countFturCovInfoPath = "$outDir/preCalInfo/countFturCovPerNt.txt" if ($countFturCovInfoPath eq "auto");
		open (TEST, "$intronBoundCovPath") || die "Cannot open intronBoundCovPath: $!\n"; close TEST;
		open (TEST, "$countFturCovInfoPath") || die "Cannot open countFturCovInfoPath: $!\n"; close TEST;
	}

	if ($unqPileupPath ne "no") {
		open (TEST, "$unqPileupPath") || die "Cannot open $unqPileupPath: $!\n"; close TEST;
	} else {
		$precalUnqLoCmPath = "$outDir/preCalInfo/precalculated.unq.txt" if ($precalUnqLoCmPath eq "auto");
		open (TEST, "$precalUnqLoCmPath") || die "Cannot open $precalUnqLoCmPath: $!\n"; close TEST;
	}

	#chop $outDir if ($outDir =~ m/\/$/); #---remove the last slash
	system "mkdir -p -m 777 $outDir";
	system "mkdir -p -m 777 $outDir/plotData";
	system "mkdir -p -m 777 $outDir/plotPdf";
	system "mkdir -p -m 777 $outDir/flankSeq";
	system "mkdir -p -m 777 $outDir/scoreReadNum";
	system "mkdir -p -m 777 $outDir/altSpliceSiteShift";
	system "mkdir -p -m 777 $outDir/preCalInfo";
	system "mkdir -p -m 777 $outDir/junctInfo";
	system "mkdir -p -m 777 $outDir/correlation";
	system "mkdir -p -m 777 $outDir/GFFAndBED";
	system "mkdir -p -m 777 $outDir/altSpliceInfo";
	system "mkdir -p -m 777 $outDir/geneBasedInfo";
	system "mkdir -p -m 777 $outDir/tmpLog";
	system "mkdir -p -m 777 $outDir/exonSkip";
	system "mkdir -p -m 777 $outDir/ovrlpAntisenseJunct";
	system "mkdir -p -m 777 $outDir/splicingEffInterval/";
	system "mkdir -p -m 777 $outDir/splicingEffInterval/fasta";
	system "mkdir -p -m 777 $outDir/splicingEffInterval/entropy";
	system "mkdir -p -m 777 $outDir/splicingEffInterval/pdf";
	system "mkdir -p -m 777 $outDir/splicingEffInterval/log";
	system "mkdir -p -m 777 $outDir/dreme/fasta/";
	
	return ($refGffPath, $canJunctionBedPath, $nonCanJunctionBedPath, $pileupPath, $minSplicingRatio, $intronBoundCovPath, $minMajorSupportRead, $geneticCodePath, $cntgFastaPath, $geneticCodeID, $ATGOnly, $unqPileupPath, $boundWidth, $missModJunctComPath, $precalUnqLoCmPath, $countFturCovInfoPath, $pctCutLim, $strndSpecificPileup, $lowComRegionPath, \%ctrgyToCountHsh, $outDir);
}
########################################################################## findFturOverlap
sub findFturOverlap {

#					The 7 scenes of overlapping and proximity 
#
#
#     case 0: complete overlapp (($refStart == $qryStart) && ($refEnd == $qryEnd))
#			
#     case 1: overlapHead    case 2: overlapTail      case 3: cover		      case 4: within		case 5: prxmtyHead	     case 6: prxmtyTail
#
#r	     |--------|		        |---------|	        |-------------|	              |-----|			    |-----|				       	                |-------|
#q	<=========>	           <==========>		          <=========>	            <==========>		     	    <==========>	      <==========>
#
#   ($refStart<$qryStart)&&	($refStart>=$qryStart)&&  ($refStart<$qryStart)&&  ($refStart>$qryStart)&&   ($refStart<$qryStart)&&	   ($refEnd>$qryStart)&&
#   ($refEnd>=$qryStart)&&	($refStart<=$qryEnd)&&	  ($refEnd>$qryEnd)	       ($refEnd<$qryEnd)	     ($refStart<$qryEnd)		   ($refStart>$qryEnd)
#   ($refEnd<=$qryEnd)	    ($refEnd>$qryEnd)												 
#

	my %refRangeXStrndHsh = %{$_[0]};
	my %qryRangeXStrndHsh = %{$_[1]}; 
	my %refStrndByFturHsh = %{$_[2]};
	my %qryStrndByFturHsh = %{$_[3]};
	my $reportExactOverlap = $_[4]; #----yes or no, if no, exactly overlap i.e. case 0 will not be reported. This is designed to prevent self-hit if reference and query are the same hash;
	my $verbose = $_[5]; #----y or n
	my $refStrndFilter = $_[6];#--- only the indicated strd will be process, use all to indicate all
	my $qryStrndFilter = $_[7];#--- only the indicated strd will be process, use all to indicate all
	my $strToPrint = $_[8];
	
	print "\n".$strToPrint."\n";
	
	my (%SSHitByRefHsh, %SSHitByQryHsh, %XStrndHitByRefHsh, %XStrndHitByQryHsh);

	foreach my $cntg (sort {$a cmp $b} keys %refRangeXStrndHsh) {
		print "Finding overlapping features on $cntg\r" if ($verbose eq "y");
		if (exists $qryRangeXStrndHsh{$cntg}) {#---if there are ftur on the $strd of $cntg of qryGffPath
			foreach my $refFtur (sort {$a cmp $b} keys %{$refRangeXStrndHsh{$cntg}}) {#--- all ftur on the $strd of $cntg of refGff
				
				#---skip if the strd of the refFtur doesnt match the filter
				next if (($refStrndFilter ne "all") and ($refStrndByFturHsh{$refFtur} ne $refStrndFilter));
				
				my $refStart = ${${$refRangeXStrndHsh{$cntg}}{$refFtur}}{"start"};
				my $refEnd = ${${$refRangeXStrndHsh{$cntg}}{$refFtur}}{"end"};
				foreach  my $qryFtur (sort {$a cmp $b} keys %{$qryRangeXStrndHsh{$cntg}}) {#--- all ftur on the $strd of $cntg of qryGffPath
					
					#---skip if the strd of the refFtur doesnt match the filter
					next if (($qryStrndFilter ne "all") and ($qryStrndByFturHsh{$qryFtur} ne $qryStrndFilter));

					my $sameStrnd = "y";
					
					#---two futrs are not on the same strd and strd has to be + or -;  (i.e. not . or * for no strand with be hit with both + or - in the SS manner)
					if (($refStrndByFturHsh{$refFtur} ne $qryStrndByFturHsh{$qryFtur}) and (($qryStrndByFturHsh{$qryFtur} eq "+") or ($qryStrndByFturHsh{$qryFtur} eq "-")) and (($refStrndByFturHsh{$refFtur} eq "+") or ($refStrndByFturHsh{$refFtur} eq "-"))) {
						$sameStrnd = "n";
					}
					
					my $qryStart = ${${$qryRangeXStrndHsh{$cntg}}{$qryFtur}}{"start"};
					my $qryEnd = ${${$qryRangeXStrndHsh{$cntg}}{$qryFtur}}{"end"};
					
					if  (($refStart == $qryStart) && ($refEnd == $qryEnd)) {#---scene 0

						if ($reportExactOverlap eq "y") {#----to prevent selfhit when ref and qry is the same hash.
							${$XStrndHitByRefHsh{$refFtur}}{$qryFtur} = 0;
							${$XStrndHitByQryHsh{$qryFtur}}{$refFtur} = 0;

							if ($sameStrnd eq "y") {	
								${$SSHitByRefHsh{$refFtur}}{$qryFtur} = 0;
								${$SSHitByQryHsh{$qryFtur}}{$refFtur} = 0;
							}

						} else {

							if ($sameStrnd eq "n") {#---no on the same strand	
								${$XStrndHitByRefHsh{$refFtur}}{$qryFtur} = 0;
								${$XStrndHitByQryHsh{$qryFtur}}{$refFtur} = 0;
							}
						}
						
					} elsif (($refStart<=$qryStart)&&($refEnd>=$qryStart)&&($refEnd<=$qryEnd)) {#---scene 1		

						${$XStrndHitByRefHsh{$refFtur}}{$qryFtur} = 1;
						${$XStrndHitByQryHsh{$qryFtur}}{$refFtur} = 1;

						if ($sameStrnd eq "y") {
							${$SSHitByRefHsh{$refFtur}}{$qryFtur} = 1;
							${$SSHitByQryHsh{$qryFtur}}{$refFtur} = 1;
						}
											
					} elsif (($refStart>=$qryStart)&&($refStart<=$qryEnd)&&($refEnd>=$qryEnd)) {#---scene 2					

						${$XStrndHitByRefHsh{$refFtur}}{$qryFtur} = 2;
						${$XStrndHitByQryHsh{$qryFtur}}{$refFtur} = 2;

						if ($sameStrnd eq "y") {
							${$SSHitByRefHsh{$refFtur}}{$qryFtur} = 2;
							${$SSHitByQryHsh{$qryFtur}}{$refFtur} = 2;
						}
					
					} elsif (($refStart<=$qryStart)&&($refEnd>=$qryEnd)) {#---scene 3		

						${$XStrndHitByRefHsh{$refFtur}}{$qryFtur} = 3;
						${$XStrndHitByQryHsh{$qryFtur}}{$refFtur} = 3;

						if ($sameStrnd eq "y") {
							${$SSHitByRefHsh{$refFtur}}{$qryFtur} = 3;
							${$SSHitByQryHsh{$qryFtur}}{$refFtur} = 3;
						}
					
					} elsif (($refStart>=$qryStart)&&($refEnd<=$qryEnd)) {#---scene 4						

						${$XStrndHitByRefHsh{$refFtur}}{$qryFtur} = 4;
						${$XStrndHitByQryHsh{$qryFtur}}{$refFtur} = 4;

						if ($sameStrnd eq "y") {
							${$SSHitByRefHsh{$refFtur}}{$qryFtur} = 4;
							${$SSHitByQryHsh{$qryFtur}}{$refFtur} = 4;
						}
					
					} elsif (($refStart<$qryStart)&&($refStart<$qryEnd)) {#---scene 5
						#---non-overlapping, nothing to be done;						
					} elsif (($refEnd>$qryStart)&&($refStart>$qryEnd)) {#---scene 6
						#---non-overlapping, nothing to be done;						
					} else {#---BUG! possibly other scene?
						die "Unexpected overlapping scene between $refFtur and $qryFtur. It's a Bug. Program qutting.\n";
					}
				}
			}
		}
	}
	
	return (\%SSHitByRefHsh, \%SSHitByQryHsh, \%XStrndHitByRefHsh, \%XStrndHitByQryHsh);
}
########################################################################## readJunctBED
sub readJunctBED {

	my $junctBEDPath = $_[0];
	
	my (%junctRngHsh, %NGSJunctStrndHsh, %junctCtngHsh, %junctReadNumHsh, %junctScoreHsh, %junctIntronRngHsh, %NGSJunctBEDCoverRngHsh);
	
	my $dupJunctNum = 0;
	
	open (INFILE, "$junctBEDPath");
	open (TMPLOG, ">$outDir/tmpLog/junctionInBedMoreThanOnce.txt");
	print "Reading $junctBEDPath\n";
	while (my $theLine = <INFILE>) {
		chomp $theLine;
		next if ($theLine =~ m/^track name/);
		my @theLineSplt = split (/\t/, $theLine);
		my $cntg = $theLineSplt[0];
		my $bedStart = $theLineSplt[1];
		my $bedEnd = $theLineSplt[2];
		my $strd = $theLineSplt[5]; #--- will be an asterick for non-cannoical
		my @blkSizesSplt = split /,/, $theLineSplt[10];
		my $blk1Size = $blkSizesSplt[0];
		my $blk2Size = $blkSizesSplt[1];
		my @readNumAndScoreSplt = split /\|/, $theLineSplt[3];
		my $readNum = substr ($readNumAndScoreSplt[0], index ($readNumAndScoreSplt[0], "=")+1);
		my $score = substr ($readNumAndScoreSplt[1], index ($readNumAndScoreSplt[1], "=")+1);
		
		my $intronStart = $bedStart + $blk1Size + 1;
		my $intronEnd = $bedEnd - $blk2Size;

		my $junctStr = $cntg.":".$intronStart.":".$intronEnd; #---assumed to be unique
		
		${${$junctIntronRngHsh{$cntg}}{$junctStr}}{"start"} = $intronStart;
		${${$junctIntronRngHsh{$cntg}}{$junctStr}}{"end"} = $intronEnd;
		
		#---multiple $junctStr may exist in HMMSplicer as unique and duplicated individual junctions are collapsed seperately
		if (not exists $junctCtngHsh{$junctStr}) { #---most of the case
			$junctCtngHsh{$junctStr} = $cntg;
			${${$junctRngHsh{$cntg}}{$junctStr}}{"start"} = $bedStart;
			${${$junctRngHsh{$cntg}}{$junctStr}}{"end"} = $bedEnd;
			$NGSJunctBEDCoverRngHsh{$junctStr} = (($intronStart - $bedStart) + ($bedEnd - $intronEnd));
			$NGSJunctStrndHsh{$junctStr} = $strd;
			$junctReadNumHsh{$junctStr} = $readNum;
			$junctScoreHsh{$junctStr} = $score;

		} else { #---appeared twice, extend the range 
			
			$dupJunctNum++;
			
			my $storedStart = ${${$junctRngHsh{$cntg}}{$junctStr}}{"start"}; 
			my $storedEnd = ${${$junctRngHsh{$cntg}}{$junctStr}}{"end"};
			my $storedScore = $junctScoreHsh{$junctStr};
			my $storedReadNum = $junctReadNumHsh{$junctStr};
			
			$junctReadNumHsh{$junctStr} = $readNum + $storedReadNum;
			$junctScoreHsh{$junctStr} = $score if ($score > $storedScore);;
			
			${${$junctRngHsh{$cntg}}{$junctStr}}{"start"} = $bedStart if ($bedStart < $storedStart);
			${${$junctRngHsh{$cntg}}{$junctStr}}{"end"} = $bedEnd if ($bedEnd > $storedEnd);
			
			$NGSJunctBEDCoverRngHsh{$junctStr} = (($intronStart - $bedStart) + ($bedEnd - $intronEnd));#---ensure the longest range

			print TMPLOG "$junctStr appeared in the bed file twice. Collpasing the two BED lines. Score = $storedScore|$score, read=$storedReadNum|$readNum\n";
			
		}

	}	
	close INFILE;
	close TMPLOG;
	
	my $junctNum = keys %NGSJunctStrndHsh;
	
	print "Totally $junctNum junctions have been stored. With $dupJunctNum of them appeared twice and collapsed.\n";

	return (\%NGSJunctStrndHsh, \%junctRngHsh, \%junctCtngHsh, \%junctScoreHsh, \%junctReadNumHsh, \%junctIntronRngHsh, \%NGSJunctBEDCoverRngHsh);
	
}
########################################################################## printCMDLogOrFinishMessage
sub printCMDLogOrFinishMessage {

	my $CMDLogOrFinishMessage = $_[0];
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $scriptNameXext = $0;
		$scriptNameXext =~ s/\.\w+$//;
		open (CMDLOG, ">>$scriptNameXext.cmd.log.txt"); #---append the CMD log file
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print CMDLOG "[".$runTime."]\t"."perl $0 ".(join " ", @ARGV)."\n";
		close CMDLOG;
		print "\n=========================================================================\n";
		print "$0 starts running at [$runTime]\n";
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print "\n=========================================================================\n";
		print "$0 finished running at [$runTime]\n";
		print "=========================================================================\n\n";
	}
	
}
########################################################################## readGff
sub readGff {

	my $refGff = $_[0];
	my %ctrgyToStoreHsh = %{$_[1]};
	
	#---the whole subrountine was inherited from pileupCounter_v0.8 so it looks a bit redundant. Will come back later to clean it up.

	#---variables to retun
	my (%nameByGeneHsh, %strndByGeneHsh, %cntgByGeneHsh, %exonRngByGeneHsh, %exonNumByCntgHsh, %geneExonLenHsh, %geneCDSLenHsh, %ctgryReadCountHsh, %CDSRngByGeneHsh, %geneByCtgryHsh, %ctgryByGeneHsh);
	my (%geneByRNAHsh, %CDSCountHsh, %exonCountHsh, %geneExonLocationHsh);
	my (%SSRngByCntgByWholeGeneHsh, %XStrndRngByCntgByWholeGeneHsh);

	#---read the gff
	open (INFILE, $refGff) || die "Cannot open $refGff";
	print "Reading $refGff for categorizing the features.\n";
	while (my $theLine = <INFILE>) {
		chomp $theLine;
		if (($theLine !~ m/^\#|^\@/) and ($theLine !~ m/\tsupercontig\t/)) {
			my @theLineSplt = split (/\t/, $theLine);
			my $seq = $theLineSplt[0];

			my $geneCategory = $theLineSplt[2];
			my $featureStart = $theLineSplt[3];
			my $featureEnd = $theLineSplt[4];
			my $geneStrnd = $theLineSplt[6];
			my $dscrptns = $theLineSplt[8];
			my @dscrptnsSplt = split /;/, $dscrptns;
			my ($unqID, $parent, $geneName);
			foreach my $theDscptn (@dscrptnsSplt) {
				if ($theDscptn =~ m/^ID=/) {$unqID = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^Parent=/) {$parent = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^description=/) {$geneName = substr ($theDscptn, index ($theDscptn, "=")+1);}
			}
	
			if ($geneCategory eq "gene") {#---gene
				
				my $geneID = $unqID;
				$geneName = $geneID if (not defined $geneName);
				$nameByGeneHsh{$geneID} = $geneName;
				$strndByGeneHsh{$geneID} = $geneStrnd;
				$cntgByGeneHsh{$geneID} = $seq;
				${${${$SSRngByCntgByWholeGeneHsh{$seq}}{$geneStrnd}}{$geneID}}{"start"} = $featureStart;
				${${${$SSRngByCntgByWholeGeneHsh{$seq}}{$geneStrnd}}{$geneID}}{"end"} = $featureEnd;

				${${$XStrndRngByCntgByWholeGeneHsh{$seq}}{$geneID}}{"start"} = $featureStart;
				${${$XStrndRngByCntgByWholeGeneHsh{$seq}}{$geneID}}{"end"} = $featureEnd;

			} elsif ($geneCategory eq "CDS") {#---Only for coding genes
				
				# The CDS is ignored at the moment, until it reaches the point that we are looking at UTRs
				#
				#my $mRNAID = $parent;
				#my $geneID = $geneByRNAHsh{$mRNAID};
				#$CDSCountHsh{$geneID}++;
				#my $CDSCount = $CDSCountHsh{$geneID};
				#${${$CDSRngByGeneHsh{$geneID}}{$CDSCount}}{"start"} = $featureStart;
				#${${$CDSRngByGeneHsh{$geneID}}{$CDSCount}}{"end"} = $featureEnd;
			 	#$geneCDSLenHsh{$geneID} = 0 if $CDSCount == 1; #---define the length hashfor the 1st time
			 	#$geneCDSLenHsh{$geneID} += ${${$CDSRngByGeneHsh{$geneID}}{$CDSCount}}{"end"} - ${${$CDSRngByGeneHsh{$geneID}}{$CDSCount}}{"start"};
				
			} elsif ($geneCategory eq "exon") {#---exon, may be exons of alternative transcripts, wiull sort out later
				my $exonID = $unqID;
				my $RNAID = $parent;
				my $geneID = $geneByRNAHsh{$RNAID};

				${${$exonRngByGeneHsh{$geneID}}{$exonID}}{"start"} = $featureStart;
				${${$exonRngByGeneHsh{$geneID}}{$exonID}}{"end"} = $featureEnd;
				
			} else {#---can be tRNA, rRNA, mRNA, repRNA, ncRNA
				my $RNAID = $unqID;
				my $geneID = $parent;
				$geneByRNAHsh{$RNAID} = $geneID;
				$ctgryByGeneHsh{$geneID} = $geneCategory;
				
				if (not(exists $ctgryReadCountHsh{$geneCategory})) {#---initialize the $geneCategory for all category 
					${$ctgryReadCountHsh{$geneCategory}}{"s"} = 0;
					${$ctgryReadCountHsh{$geneCategory}}{"a"} = 0;
				}
			}
		}#---end of if (($theLine !~ m/^\#|^\@/) and ($theLine !~ m/\tsupercontig\t/)) {
	}#---end of while (my $theLine = <INFILE>)
	close INFILE;
	
	#---remove the unwanted ctgry
	if (not exists $ctrgyToStoreHsh{"all"}) {
		foreach my $geneID (keys %strndByGeneHsh) {
			my $ctgry = $ctgryByGeneHsh{$geneID};
			my $strnd = $strndByGeneHsh{$geneID};
			my $cntg = $cntgByGeneHsh{$geneID};
			if (not exists $ctrgyToStoreHsh{$ctgry}) {
				delete $ctgryByGeneHsh{$geneID};
				delete $strndByGeneHsh{$geneID};
				delete $cntgByGeneHsh{$geneID};
				delete $exonRngByGeneHsh{$geneID};
				delete ${${$SSRngByCntgByWholeGeneHsh{$cntg}}{$strnd}}{$geneID};
				delete ${$XStrndRngByCntgByWholeGeneHsh{$cntg}}{$geneID};
				delete $nameByGeneHsh{$geneID};
			}
		}
	}

	my ($refIntronRngByGeneHsh_ref, $refJunctStrByIntronIDHsh_ref, $refIntronRngXS_ref, $refStrndByIntronHsh_ref, $geneIDByJunctStrHsh_ref) = getIntronFromExonRng(\%exonRngByGeneHsh, \%cntgByGeneHsh, \%strndByGeneHsh, "Getting the intron boundaries on refGff");
	
	return (\%SSRngByCntgByWholeGeneHsh, \%ctgryByGeneHsh, \%strndByGeneHsh, \%XStrndRngByCntgByWholeGeneHsh, \%exonRngByGeneHsh, $refIntronRngByGeneHsh_ref, $refJunctStrByIntronIDHsh_ref, $refIntronRngXS_ref, $refStrndByIntronHsh_ref, \%nameByGeneHsh, $geneIDByJunctStrHsh_ref, \%cntgByGeneHsh);
}
########################################################################## defineOverlappingJunctBEDRng
sub defineOverlappingJunctBEDRng {

	my %SSOvrlpJunctBEDRngHsh = %{$_[0]};
	my %junctRngHsh = %{$_[1]};
	my %NGSJunctCntgHsh = %{$_[2]};
	my %NGSJunctStrndHsh = %{$_[3]};
	
	#--- define the clusters
	my ($junctBEDRngClusterNum, $junctBEDRngClusterNameByJunctHsh_ref, $junctBEDRngClusterJunctHsh_ref) = recursivelyFindJunctionCluster(\%SSOvrlpJunctBEDRngHsh, "junction BED range", "bc_");
	my %junctBEDRngClusterJunctHsh = %{$junctBEDRngClusterJunctHsh_ref};
	my %junctBEDRngClusterNameByJunctHsh = %{$junctBEDRngClusterNameByJunctHsh_ref};

	#--- assign a cluster to the orphan
	my $orphanBEDRngNum = 0;
	
	foreach my $junctStr (sort {$a cmp $b} keys %NGSJunctStrndHsh) {
		$orphanBEDRngNum++;
		my $clusterName = "bo_".$orphanBEDRngNum;
		if (not exists $junctBEDRngClusterNameByJunctHsh{$junctStr}) { #--it is an orphan
			$junctBEDRngClusterNameByJunctHsh{$junctStr} = $clusterName;
			push @{$junctBEDRngClusterJunctHsh{$clusterName}}, $junctStr;
		}
	}
	
	print "Defining the junction BED range clusters\n";
	#---define the cluster range
	my (%junctBEDClusterRngHsh, %junctBEDClusterStrndHsh, %junctBEDClusterCntgHsh);
	foreach my $clusterName (sort {$a cmp $b} keys %junctBEDRngClusterJunctHsh) {
		my ($cntg, $strd, @tmpBoundAry);
		foreach my $junctStr (@{$junctBEDRngClusterJunctHsh{$clusterName}})	{
			$cntg = $NGSJunctCntgHsh{$junctStr};
			$strd = $NGSJunctStrndHsh{$junctStr};
			push @tmpBoundAry, ${${$junctRngHsh{$cntg}}{$junctStr}}{"start"};
			push @tmpBoundAry, ${${$junctRngHsh{$cntg}}{$junctStr}}{"end"};
		}
		my @tmpBoundSortedAry = sort {$a <=> $b} @tmpBoundAry;
		
		${${$junctBEDClusterRngHsh{$cntg}}{$clusterName}}{"start"} = $tmpBoundSortedAry[0];
		${${$junctBEDClusterRngHsh{$cntg}}{$clusterName}}{"end"} = $tmpBoundSortedAry[-1];
		$junctBEDClusterStrndHsh{$clusterName} = $strd;
		$junctBEDClusterCntgHsh{$clusterName} = $cntg;
		#print $clusterName."\t".${${$junctBEDClusterRngHsh{$cntg}}{$clusterName}}{"start"}."\t".${${$junctBEDClusterRngHsh{$cntg}}{$clusterName}}{"end"}."\t".$strd."\t".$cntg."\n";
	}
	
	return (\%junctBEDRngClusterJunctHsh, \%junctBEDRngClusterNameByJunctHsh, \%junctBEDClusterRngHsh, \%junctBEDClusterStrndHsh, \%junctBEDClusterCntgHsh);
}
########################################################################## recursivelyFindJunctionCluster
sub recursivelyFindJunctionCluster {
	
	#---works for both intron range and BED range
	
	our %SSOvrlpRngHsh = %{$_[0]};
	my $reportString = $_[1];
	my $clusterIDTag = $_[2];
	
	#---var to return
	our %clusterJunctRngHsh = ();
	our %jClusterNameHsh = (); 

	#---check each transfrag
	my $i = my $clusterNum = 0;
	
	print "Start clustering the $reportString.\n";
	
	foreach my $refJunct (keys %SSOvrlpRngHsh) {
	
		if (not exists $jClusterNameHsh{$refJunct}) {
			$clusterNum++; $i++;
			my $clusterName = $clusterIDTag.$clusterNum;
			my @intitialJunctRecursiveAry = ($refJunct);
			my $jClusterNameHsh_ref = \%jClusterNameHsh;
			my $clusterJunctRngHsh_ref = \%clusterJunctRngHsh;
	
			recursiveClusterExtension(\@intitialJunctRecursiveAry, $clusterName);
		
			if ($i eq 1000) {
				$i = 0;	
				print $clusterNum." $reportString clusters identified.\n";
			}
		}
	}

	return ($clusterNum, \%jClusterNameHsh, \%clusterJunctRngHsh);
	
	####################################################recursiveClusterExtension####################################################
	sub recursiveClusterExtension {
	
		my @recusiveJunctInAry = @{$_[0]};
		my $r_clusterName = $_[1];

		my @recusiveJunctOutAry;
		
		foreach my $recusiveJunctIn (@recusiveJunctInAry) {

			if (not exists $jClusterNameHsh{$recusiveJunctIn}) { 
				$jClusterNameHsh{$recusiveJunctIn} = $r_clusterName;
				push @{$clusterJunctRngHsh{$r_clusterName}}, $recusiveJunctIn;

				foreach my $r_refJunct (keys %{$SSOvrlpRngHsh{$recusiveJunctIn}}) {#---check each 
					
					if (not exists $jClusterNameHsh{$r_refJunct}) { 
						$jClusterNameHsh{$r_refJunct} = $r_clusterName;
						push @{$clusterJunctRngHsh{$r_clusterName}}, $r_refJunct;
				
						#----look for the junctions hit this junction
						foreach my $r_qryJunct (keys %{$SSOvrlpRngHsh{$r_refJunct}}) {
					
							#---only if the hit is recusiveJunctIn and the hit has not been reached before
							if (($r_qryJunct ne $recusiveJunctIn) and (not exists $jClusterNameHsh{$r_qryJunct}) and ($r_qryJunct ne $r_refJunct)) {
								push @recusiveJunctOutAry, $r_qryJunct;
							}
						}
					}
				}
			}
		}

		#---if there are junct hit, call itself again
		recursiveClusterExtension(\@recusiveJunctOutAry, $r_clusterName) if (@recusiveJunctOutAry > 0);
	}	
}
########################################################################## scanJunctionIntronBoundCoverage
sub scanJunctionIntronBoundCoverage {#----scan for the coverage flanking the intron bounds and record the splicing ratio 

	my %junctIntronRngHsh = %{$_[0]};
	my %fturCountRngXSHsh = %{$_[1]};
	my %fturCountExonRngByGeneHsh = %{$_[2]};

	#---var to return
	my %junctIntronBoundCovHsh = my %countFturCovInfoHsh = ();

	#---read or calculate the intron bound coverage
	if (($intronBoundCovPath ne "no") and ($countFturCovInfoPath ne "no")){#---read
		
		print "Reading $intronBoundCovPath for intron bound coverage\n";
		
		open (INTRNBOUNDCOV, "$intronBoundCovPath")|| die "Can't read $intronBoundCovPath :$!\n";
		
		while (my $theLine = <INTRNBOUNDCOV>) {
			chomp $theLine;
			my @theLineSplt = split /\t/, $theLine;
			my $junctStr = $theLineSplt[0];
			${${$junctIntronBoundCovHsh{$junctStr}}{"start"}}{"+"} = $theLineSplt[1];
			${${$junctIntronBoundCovHsh{$junctStr}}{"start"}}{"-"} = $theLineSplt[2];
			${${$junctIntronBoundCovHsh{$junctStr}}{"end"}}{"+"} = $theLineSplt[3];
			${${$junctIntronBoundCovHsh{$junctStr}}{"end"}}{"-"} = $theLineSplt[4];
		}
		close INTRNBOUNDCOV;

		print "Reading $countFturCovInfoPath for coverage info of countFtur\n";
		
		open (CNTFTURCOVINFO, "$countFturCovInfoPath")|| die "Can't read $countFturCovInfoPath :$!\n";
		
		while (my $theLine = <CNTFTURCOVINFO>) {
			chomp $theLine;
			my @theLineSplt = split /\t/, $theLine;
			my $geneID = $theLineSplt[0];
			${$countFturCovInfoHsh{$geneID}}{"length"} = $theLineSplt[1];
			${$countFturCovInfoHsh{$geneID}}{"plusCov"} = $theLineSplt[2]; 
			${$countFturCovInfoHsh{$geneID}}{"minusCov"} = $theLineSplt[3];
			${$countFturCovInfoHsh{$geneID}}{"bothCov"} = $theLineSplt[2] + $theLineSplt[3];
		}
		close CNTFTURCOVINFO;
		
	} else {#calculate
		
		open (PILEUPFILE, "$pileupPath");
		open (INTRNBOUNDCOV, ">$outDir/preCalInfo/intronBoundCov.txt");
		open (CNTFTURCOVINFO, ">$outDir/preCalInfo/countFturCovPerNt.txt");
	
		print "Reading $pileupPath.\n";
	
		my %tmpCovByPosHsh = ();
		my $flankSize = 30;
		
		#--get the first line
		my $theCurntLine;
		while ($theCurntLine = <PILEUPFILE>) {
			if ($theCurntLine !~ m/^[\@|\#]/) {#---ignore comments
				chomp $theCurntLine;
				last;
			}
		}
		
		#--go through the rest of the file
		while (my $theNextLine = <PILEUPFILE>) {
			next if ($theNextLine =~ m/^[\@|\#]/);
			chomp $theNextLine;
		
			my @theNextLineSplt = split /\t/, $theNextLine;
			my $nextCtng = $theNextLineSplt[0];
	
			my @theCurntLineSplt = split /\t/, $theCurntLine;
			my $curntCtng = $theCurntLineSplt[0];
	
			my $pos = $theCurntLineSplt[1];
			my $plusCov = $theCurntLineSplt[3];
			my $minusCov = $theCurntLineSplt[4];
		
			#---store the cov if plus or minus strand > 0;
			if (($plusCov > 0) or ($minusCov > 0)) {
				${$tmpCovByPosHsh{$pos}}{"+"} = $plusCov;
				${$tmpCovByPosHsh{$pos}}{"-"} = $minusCov;
			}
			
			if (eof(PILEUPFILE)) {#---last line of the file
				my $pos = $theNextLineSplt[1];
				my $plusCov = $theNextLineSplt[3];
				my $minusCov = $theNextLineSplt[4];
	
				#---store the cov if plus or minus strand > 0;
				if (($plusCov > 0) or ($minusCov > 0)) {
					${$tmpCovByPosHsh{$pos}}{"+"} = $plusCov;
					${$tmpCovByPosHsh{$pos}}{"-"} = $minusCov;
				}
			}
			
			#---change contig or end of file
			if (($curntCtng ne $nextCtng) or (eof(PILEUPFILE))) {
				
				print "Scanning intron bound coverage and refGene in contig $curntCtng.\n";
				
				foreach my $geneID (keys %{$fturCountRngXSHsh{$curntCtng}}) {
					my $plusCovSum = my $minusCovSum = my $totalPos = 0;

					foreach my $exonID (keys %{$fturCountExonRngByGeneHsh{$geneID}}) {
						
						my $exonStart = ${${$fturCountExonRngByGeneHsh{$geneID}}{$exonID}}{"start"};
						my $exonEnd = ${${$fturCountExonRngByGeneHsh{$geneID}}{$exonID}}{"end"};
						
						for my $pos ($exonStart..$exonEnd) {
							$totalPos++;
							if (exists $tmpCovByPosHsh{$pos}) {
								$plusCovSum += ${$tmpCovByPosHsh{$pos}}{"+"};
								$minusCovSum += ${$tmpCovByPosHsh{$pos}}{"-"};
							}
						}
					}
					
					${$countFturCovInfoHsh{$geneID}}{"length"} = $totalPos;
					${$countFturCovInfoHsh{$geneID}}{"plusCov"} = $plusCovSum; 
					${$countFturCovInfoHsh{$geneID}}{"minusCov"} = $minusCovSum;
					
					print CNTFTURCOVINFO $geneID."\t".$totalPos."\t".$plusCovSum."\t".$minusCovSum."\n";
				}
				
				foreach my $junctStr (keys %{$junctIntronRngHsh{$curntCtng}}) {
					
					my $intronStart = ${${$junctIntronRngHsh{$curntCtng}}{$junctStr}}{"start"};
					my $intronEnd = ${${$junctIntronRngHsh{$curntCtng}}{$junctStr}}{"end"};
					
					my $startPlusCovSum = my $startMinusCovSum = my $endPlusCovSum = my $endMinusCovSum = 0;
					
					#---intron start position flank region cov
					for my $pos (($intronStart)..($intronStart + $flankSize)) {
						if (exists $tmpCovByPosHsh{$pos}) {
							$startPlusCovSum += ${$tmpCovByPosHsh{$pos}}{"+"};
							$startMinusCovSum += ${$tmpCovByPosHsh{$pos}}{"-"};
						}
					}
					
					${${$junctIntronBoundCovHsh{$junctStr}}{"start"}}{"+"} = $startPlusCovSum/$flankSize;
					${${$junctIntronBoundCovHsh{$junctStr}}{"start"}}{"-"} = $startMinusCovSum/$flankSize;
	
					#---intron end position flank region cov
					for my $pos (($intronEnd - $flankSize)..$intronEnd) {
						if (exists $tmpCovByPosHsh{$pos}) {
							$endPlusCovSum += ${$tmpCovByPosHsh{$pos}}{"+"};
							$endMinusCovSum += ${$tmpCovByPosHsh{$pos}}{"-"};
						}
					}
					
					${${$junctIntronBoundCovHsh{$junctStr}}{"end"}}{"+"} = $endPlusCovSum/$flankSize;
					${${$junctIntronBoundCovHsh{$junctStr}}{"end"}}{"-"} = $endMinusCovSum/$flankSize;
					
					print INTRNBOUNDCOV $junctStr."\t";
					print INTRNBOUNDCOV ${${$junctIntronBoundCovHsh{$junctStr}}{"start"}}{"+"}."\t";
					print INTRNBOUNDCOV ${${$junctIntronBoundCovHsh{$junctStr}}{"start"}}{"-"}."\t";
					print INTRNBOUNDCOV ${${$junctIntronBoundCovHsh{$junctStr}}{"end"}}{"+"}."\t";
					print INTRNBOUNDCOV ${${$junctIntronBoundCovHsh{$junctStr}}{"end"}}{"-"}."\n";
				} #---end of foreach my $junctStr (keys %{$junctIntronRngHsh{$curntCtng}})
				
				foreach my $junctStr (keys %{$junctIntronRngHsh{$curntCtng}}) {
				
				}
				
				
				%tmpCovByPosHsh = (); #---empty the hash
	
			} #---end of if (($curntCtng ne $nextCtng) or (eof(PILEUPFILE))) {
			
			$theCurntLine = $theNextLine;
			
		}#---end of while (my $theNextLine = <PILEUPFILE>) {
	
		close PILEUPFILE;
		close INTRNBOUNDCOV;
		close CNTFTURCOVINFO;
	}
	
	print "\n";
	
	return (\%junctIntronBoundCovHsh, \%countFturCovInfoHsh);
	
}
########################################################################## defineIntronClusters
sub defineIntronClusters {
	
	our %SSOvrlpJunctIntronHsh = %{$_[0]};
	our %NGSJunctStrndHsh = %{$_[1]};
	our %junctScoreHsh = %{$_[2]};
	our %junctReadNumHsh = %{$_[3]};
	our %junctIntronBoundCovHsh = %{$_[4]};
	our %junctIntronRngHsh = %{$_[5]};
	our %XSOvrlpJunctHsh = %{$_[6]};
	our %refGeneIDByJunctStrHsh = %{$_[7]};
	our %NGSJunctUnqHsh = %{$_[8]};
	our %NGSJunctSeqHsh = %{$_[9]};
	our %NGSJunctBEDCoverRngHsh = %{$_[10]};
	
	#---define splcing raio and the major junctions
	my ($junctSplicingRatioHsh_ref, $majorJunctAvgSplcingRatioHsh_ref, $allJunctAvgSplcingRatioHsh_ref) = defineMajorJunctions();

	#---define superJunctions
	my $superJunctHsh_ref = defineSuperJunctions();
	
	#---define overlapping intron clusters excluding superJunctions
	my ($clusterAllJunctHsh_ref, $jClusterNameByJunctHsh_ref, $jClusterInfoHsh_ref, $jClusterSSOverlapHsh_ref, $NGSJunctInfoHsh_ref) = defineIntronCluster($superJunctHsh_ref, $majorJunctAvgSplcingRatioHsh_ref, $allJunctAvgSplcingRatioHsh_ref);
	
	#---define the exon skipping events
	my ($totalExonSkippingTypeHsh_ref, $majorExonSkippingTypeHsh_ref, $exactExonSkippingClusterHsh_ref, $superJunctOvrlpClusterHsh_ref) = defineExonSkipping($superJunctHsh_ref, $majorJunctAvgSplcingRatioHsh_ref, $jClusterNameByJunctHsh_ref);
	
	#---define the alternative 5' and 3' splicing sites
	my ($splicingSiteDiffAllHsh_ref, $splicingSiteDiffHshByClusterHsh_ref);
	($splicingSiteDiffAllHsh_ref, $splicingSiteDiffHshByClusterHsh_ref, $NGSJunctInfoHsh_ref) = defineAltSplicingSiteWithinCluster($clusterAllJunctHsh_ref, $majorJunctAvgSplcingRatioHsh_ref, $NGSJunctInfoHsh_ref);
	my %NGSJunctInfoHsh = %{$NGSJunctInfoHsh_ref};
	$NGSJunctInfoHsh_ref = \%NGSJunctInfoHsh;
	
	#---define the junctions on the prominent isoform
	my $junctOnProminentIsofmHsh_ref; 
	($jClusterInfoHsh_ref, $junctOnProminentIsofmHsh_ref) = defineJunctionsOfProminentIsoform($jClusterNameByJunctHsh_ref, $jClusterInfoHsh_ref, $jClusterSSOverlapHsh_ref);
	my %jClusterInfoHsh = %{$jClusterInfoHsh_ref}; $jClusterInfoHsh_ref = \%jClusterInfoHsh;
	
	return ($junctSplicingRatioHsh_ref, $majorJunctAvgSplcingRatioHsh_ref, $superJunctHsh_ref, $clusterAllJunctHsh_ref, $jClusterNameByJunctHsh_ref, $jClusterInfoHsh_ref, $jClusterSSOverlapHsh_ref, $totalExonSkippingTypeHsh_ref, $majorExonSkippingTypeHsh_ref, $exactExonSkippingClusterHsh_ref, $splicingSiteDiffAllHsh_ref, $splicingSiteDiffHshByClusterHsh_ref, $junctOnProminentIsofmHsh_ref, $superJunctOvrlpClusterHsh_ref, $NGSJunctInfoHsh_ref);
	
	################################################## defineMajorJunctions  #########################################################################
	sub defineMajorJunctions {
	
		print "Defining major junctions.\n";

		#---define the major junctions;
		my $histogramCutOff = 20;
		my (@endSplicingRatioToPlotAry, @startSplicingRatioToPlotAry, %startSplicingRatioToPlotHsh, %endSplicingRatioToPlotHsh, %junctSplicingRatioHsh, %majorJunctAvgSplcingRatioHsh, %allJunctAvgSplcingRatioHsh);
		
		open (TMPLOG, ">$outDir/junctInfo/majorMinorJunctBasedOnSplicingRatio.txt");
		foreach my $junctStr (keys %junctIntronBoundCovHsh) {
			my $supportReadNum = $junctReadNumHsh{$junctStr};		
			my $startBothSum =  ${${$junctIntronBoundCovHsh{$junctStr}}{"start"}}{"+"} + ${${$junctIntronBoundCovHsh{$junctStr}}{"start"}}{"-"};
			my $endBothSum = ${${$junctIntronBoundCovHsh{$junctStr}}{"end"}}{"+"} + ${${$junctIntronBoundCovHsh{$junctStr}}{"end"}}{"-"};
			my $startSplicingRatio = my $endSplicingRatio = 99999;
			$startSplicingRatio = $supportReadNum/$startBothSum if ($startBothSum > 0); 
			$endSplicingRatio = $supportReadNum/$endBothSum if ($endBothSum > 0);
			${$junctSplicingRatioHsh{$junctStr}}{"start"} = $startSplicingRatio;
			${$junctSplicingRatioHsh{$junctStr}}{"end"} = $endSplicingRatio;
	
			my $startSplicingRatioToPlot = $startSplicingRatio;
			my $endSplicingRatioToPlot = $endSplicingRatio; 
	
			$startSplicingRatioToPlot = $histogramCutOff if ($startSplicingRatioToPlot > $histogramCutOff);
			$endSplicingRatioToPlot = $histogramCutOff if ($endSplicingRatioToPlot > $histogramCutOff);
	
			$startSplicingRatioToPlotHsh{$startSplicingRatioToPlot}++;
			$endSplicingRatioToPlotHsh{$endSplicingRatioToPlot}++;
	
			push @endSplicingRatioToPlotAry, $endSplicingRatioToPlot;
			push @startSplicingRatioToPlotAry, $startSplicingRatioToPlot;
			
			my $majorJunct = "n";
			
			$allJunctAvgSplcingRatioHsh{$junctStr} = ($endSplicingRatio + $startSplicingRatio)/2;
			
			if (($startSplicingRatio >= $minSplicingRatio) and ($endSplicingRatio >= $minSplicingRatio) and ($supportReadNum >= $minMajorSupportRead)) {
				$majorJunct = "y";
				$majorJunctAvgSplcingRatioHsh{$junctStr} = sprintf "%0.2f", ($endSplicingRatio + $startSplicingRatio)/2;
			}
			print TMPLOG $junctStr."\t".$startBothSum."\t".$endBothSum."\t".$startSplicingRatio."\t".$endSplicingRatio."\t".$majorJunct."\n";
		}
		close TMPLOG;
		
		my $totalJunctNum = keys %junctIntronBoundCovHsh;
		my $majorJunctNum = keys %majorJunctAvgSplcingRatioHsh;
		
		print $majorJunctNum." out of $totalJunctNum junctions having boundary splicing ratio >= $minSplicingRatio and supporting read >= $minMajorSupportRead.\n";
		
		#---plot the splicing ratio info
		my $cumulStartSplicingRatioToPlotHsh_ref = generateCumulativePctHash(\%startSplicingRatioToPlotHsh, "a");
		my $cumulEndSplicingRatioToPlotHsh_ref = generateCumulativePctHash(\%endSplicingRatioToPlotHsh, "a");
		GNUPLOTAryHistogram(\@startSplicingRatioToPlotAry, "$outDir/plotPdf/startSplicingRatioCount.pdf", "$outDir/plotData/startSplicingRatioCount.dat", "splicing ratio", "frequency", "splicing ratio of intron start bound", "n");
		GNUPLOTAryHistogram(\@endSplicingRatioToPlotAry, "$outDir/plotPdf/endSplicingRatioCount.pdf", "$outDir/plotData/endSplicingRatioCount.dat", "splicing ratio", "frequency", "splicing ratio of intron end bound", "n");
		GNUPlotXYScatterWithLines($cumulEndSplicingRatioToPlotHsh_ref, "$outDir/plotPdf/cumulEndSplicingRatioCount.pdf", "$outDir/plotData/cumulEndSplicingRatioCount.dat", "smaller than splicing ratio", "percentage", "linear", "linear", "splicing ratio of intron end bound", "n");
		GNUPlotXYScatterWithLines($cumulStartSplicingRatioToPlotHsh_ref, "$outDir/plotPdf/cumulStartSplicingRatioCount.pdf", "$outDir/plotData/cumulStartSplicingRatioCount.dat", "smaller than splicing ratio", "percentage", "linear", "linear", "splicing ratio of intron start bound", "n");
		
		return (\%junctSplicingRatioHsh, \%majorJunctAvgSplcingRatioHsh, \%allJunctAvgSplcingRatioHsh);
	}
	################################################## defineSuperJunctions  #########################################################################
	sub defineSuperJunctions {#---A superjunction is defined as a junctiuon that is overlapping with two or more non-overlapping junctions
		
		print "Defining superJunctions.\n";

		my %superJunctHsh;
		
		open (TMPLOG, ">$outDir/junctInfo/allOverlappingJunctions.txt");
		foreach my $refJunctStr (sort {$a cmp $b} keys %SSOvrlpJunctIntronHsh) {
			print TMPLOG $refJunctStr;
			foreach my $qryJunctStr (sort {$a cmp $b} keys %{$SSOvrlpJunctIntronHsh{$refJunctStr}}) {
				print TMPLOG "\t".$qryJunctStr;
				foreach my $refOverlapJunct (keys %{$SSOvrlpJunctIntronHsh{$refJunctStr}}) {
					$superJunctHsh{$refJunctStr}++ if ((not exists ${$SSOvrlpJunctIntronHsh{$refOverlapJunct}}{$qryJunctStr}) and ($qryJunctStr ne $refOverlapJunct));
				}
			}
			print TMPLOG "\n";
		}
		close TMPLOG;
		
		print "Defining overlapping superJunctions\n";
		#----construct an hash to contain overlapping superJunctions
		my %SSOvrlpingSuperJunctHsh;
		foreach my $superJunct (sort {$a cmp $b} keys %superJunctHsh) {
			foreach my $qryJunctStr (sort {$a cmp $b} keys %{$SSOvrlpJunctIntronHsh{$superJunct}}) {
				if (exists $superJunctHsh{$qryJunctStr}) {
					${$SSOvrlpingSuperJunctHsh{$superJunct}}{$qryJunctStr} = ${$SSOvrlpJunctIntronHsh{$superJunct}}{$qryJunctStr};
				}
			}
		}

		open (TMPLOG, ">$outDir/junctInfo/overlappingSuperJunctions.txt");
		foreach my $refJunctStr (sort {$a cmp $b} keys %SSOvrlpingSuperJunctHsh) {
			print TMPLOG $refJunctStr;
			foreach my $qryJunctStr (sort {$a cmp $b} keys %{$SSOvrlpingSuperJunctHsh{$refJunctStr}}) {
				print TMPLOG "\t".$qryJunctStr;
			}
			print TMPLOG "\n";
		}
		close TMPLOG;
		
		open (TMPLOG, ">$outDir/junctInfo/allSuperJunctionsOverlappings.txt");
		foreach my $superJunct (sort {$a cmp $b} keys %superJunctHsh) {
			print TMPLOG $superJunct;
			foreach my $qryJunctStr (sort {$a cmp $b} keys %{$SSOvrlpJunctIntronHsh{$superJunct}}) {
				print TMPLOG "\t".$qryJunctStr;
			}
			print TMPLOG "\n";
		}
		close TMPLOG;
		
		return \%superJunctHsh;
	}
	############################################## defineIntronCluster ########################################################
	sub defineIntronCluster {
	
		print "Defining super-/inferior- junction clusters.\n";
	
		my %superJunctHsh = %{$_[0]};
		my %majorJunctAvgSplcingRatioHsh = %{$_[1]};
		my %allJunctAvgSplcingRatioHsh = %{$_[2]};

		#----construct an hash to contain overlapping junctions without superJunctions
		my (%SSOvrlpingJunctXSuperHsh, %SSOvrlpingJunctOnlySuperHsh);
		foreach my $refJunctStr (sort {$a cmp $b} keys %SSOvrlpJunctIntronHsh) {
			foreach my $qryJunctStr (sort {$a cmp $b} keys %{$SSOvrlpJunctIntronHsh{$refJunctStr}}) {
				if ((not exists $superJunctHsh{$refJunctStr}) and (not exists $superJunctHsh{$qryJunctStr})) {
					#---non-superJunction overlpping hash
					${$SSOvrlpingJunctXSuperHsh{$refJunctStr}}{$qryJunctStr} = ${$SSOvrlpJunctIntronHsh{$refJunctStr}}{$qryJunctStr};

				} elsif ((exists $superJunctHsh{$refJunctStr}) and (exists $superJunctHsh{$qryJunctStr})) {
					#---Only superJunction overlpping hash
					${$SSOvrlpingJunctOnlySuperHsh{$refJunctStr}}{$qryJunctStr} = ${$SSOvrlpJunctIntronHsh{$refJunctStr}}{$qryJunctStr};
				}
			}
		}
		
		my (%clusterAllJunctHsh, %jClusterNameByJunctHsh);
		
		#---recursively find the inferior clusters
		my ($inferiorClusterNum, $inferiorClusterNameByJunctHsh, $inferiorClusterJunctHsh_ref) = recursivelyFindJunctionCluster(\%SSOvrlpingJunctXSuperHsh, "inferior junction intron", "ic_");
		my %inferiorClusterJunctHsh = %{$inferiorClusterJunctHsh_ref};

		print "$inferiorClusterNum inferior junction clusters were identified.\n";

		#----transfer the inferior results to the big hash
		foreach my $clusterName (keys %inferiorClusterJunctHsh) {
			foreach my $junctStr (@{$inferiorClusterJunctHsh{$clusterName}}) {
				$jClusterNameByJunctHsh{$junctStr}  = $clusterName;
				push @{$clusterAllJunctHsh{$clusterName}}, $junctStr;
			}
			delete $inferiorClusterJunctHsh{$clusterName};
		}
		
		#---recursively find the super clusters
		my ($superClusterNum, $superClusterNameByJunctHsh, $superClusterJunctHsh_ref) = recursivelyFindJunctionCluster(\%SSOvrlpingJunctOnlySuperHsh, "super junction intron", "sc_");
		
		my %superClusterJunctHsh = %{$superClusterJunctHsh_ref};

		my %splittedSuperClusterJunctHsh;
		print "$superClusterNum super junction clusters were identified.\n";
		
		#---split superjunction cluster based on the inferior j clusters they overlap
		foreach my $clusterName (keys %superClusterJunctHsh) {

			my (%tmpInferiorJClusterHsh);
			foreach my $superJunctStr (@{$superClusterJunctHsh{$clusterName}}) {
				my %tmpJClusterHsh;
				#----record for jClusters each junctStr overlaps
				foreach my $inferiorJunctStr (keys %{$SSOvrlpJunctIntronHsh{$superJunctStr}}) {
					if (not exists $superJunctHsh{$inferiorJunctStr}) {#---inferior only, no super junct
						if (exists $jClusterNameByJunctHsh{$inferiorJunctStr}) {#---cluster
							my $inferoirJCluster = $jClusterNameByJunctHsh{$inferiorJunctStr};
							push @{$tmpInferiorJClusterHsh{$superJunctStr}}, $inferoirJCluster if (not exists $tmpJClusterHsh{$inferoirJCluster}); 
							$tmpJClusterHsh{$inferoirJCluster}++;
						} else {#---orphan
							push @{$tmpInferiorJClusterHsh{$superJunctStr}}, $inferiorJunctStr if (not exists $tmpJClusterHsh{$inferiorJunctStr}); 
							$tmpJClusterHsh{$inferiorJunctStr}++;
						}
					}
				}
			}
			#----check for different jClusters overlapping
			my %superJStrByInferiorJClusterStrHsh;
			foreach my $superJunctStr (keys %tmpInferiorJClusterHsh) {
				my @sortedJClusterAry = @{$tmpInferiorJClusterHsh{$superJunctStr}};
				@sortedJClusterAry = sort {$a cmp $b} @sortedJClusterAry;
				my $inferiorJClusterStr = join "", @sortedJClusterAry;
				push @{$superJStrByInferiorJClusterStrHsh{$inferiorJClusterStr}}, $superJunctStr;
			}
			my $subClusterNum = 0;
			foreach my $inferiorJClusterStr (keys %superJStrByInferiorJClusterStrHsh) {
				@{$splittedSuperClusterJunctHsh{$clusterName.".".$subClusterNum}} = @{$superJStrByInferiorJClusterStrHsh{$inferiorJClusterStr}};
				$subClusterNum++;
			}
		}

		#----transfer the superior results to the big hash
		foreach my $clusterName (keys %splittedSuperClusterJunctHsh) {
			foreach my $junctStr (@{$splittedSuperClusterJunctHsh{$clusterName}}) {
				$jClusterNameByJunctHsh{$junctStr}  = $clusterName;
				push @{$clusterAllJunctHsh{$clusterName}}, $junctStr;
			}
			delete $splittedSuperClusterJunctHsh{$clusterName};
		}

		#---assigne numbers for orphans also
		my $inferiorOphranNum = my $superOphranNum = 0;
		foreach my $junctStr (sort {$a cmp $b} keys %NGSJunctStrndHsh) {
			my $strand = $NGSJunctStrndHsh{$junctStr};
			if (not exists $jClusterNameByJunctHsh{$junctStr}) {#---not in the clusters, no overlapping
				if (not exists $superJunctHsh{$junctStr}) {#---non superJunction
					$inferiorOphranNum++;
					my $clusterName = "io_".$inferiorOphranNum;
					push @{$clusterAllJunctHsh{$clusterName}}, $junctStr;
					$jClusterNameByJunctHsh{$junctStr}  = $clusterName;
				} else {
					$superOphranNum++;
					my $clusterName = "so_".$superOphranNum;
					push @{$clusterAllJunctHsh{$clusterName}}, $junctStr;
					$jClusterNameByJunctHsh{$junctStr}  = $clusterName;
				}
			}
		}

		print "$superOphranNum super junction orphans were identified.\n";
		print "$inferiorOphranNum inferior junction orphans were identified.\n";

		#---get all info of the clusters
		my %jClusterInfoHsh;
		open (JCLUSBED, ">$outDir/GFFAndBED/junction.Cluster.bed");
		print JCLUSBED "track name=jCluster description=jCluster useScore=1\n";

		print "Getting all cluster information.\n";
		foreach my $clusterName (sort {$a cmp $b} keys %clusterAllJunctHsh) {
			${$jClusterInfoHsh{$clusterName}}{"majorJunctNum"} = 0;
			my $totalJNum = @{$clusterAllJunctHsh{$clusterName}};
			${$jClusterInfoHsh{$clusterName}}{"totalJunctNum"} = $totalJNum;
			my $prominentJunctStr;
			my $prominentJunctSupportReadNum = 0;
			${$jClusterInfoHsh{$clusterName}}{"majorCluster"} = "n";
			foreach my $junctStr (@{$clusterAllJunctHsh{$clusterName}}) {
				${$jClusterInfoHsh{$clusterName}}{"majorJunctNum"}++ if (exists $majorJunctAvgSplcingRatioHsh{$junctStr});			
				if ($junctReadNumHsh{$junctStr} > $prominentJunctSupportReadNum) { 
					$prominentJunctStr = $junctStr;
					$prominentJunctSupportReadNum = $junctReadNumHsh{$junctStr};
				}
			}
			
			${$jClusterInfoHsh{$clusterName}}{"majorCluster"} = "y" if (exists $majorJunctAvgSplcingRatioHsh{$prominentJunctStr});
			${$jClusterInfoHsh{$clusterName}}{"prominentJunct"} = $prominentJunctStr;
			${$jClusterInfoHsh{$clusterName}}{"strand"} = $NGSJunctStrndHsh{$prominentJunctStr};
			${$jClusterInfoHsh{$clusterName}}{"mostSupportReadNum"} = $prominentJunctSupportReadNum;
			
			#---print the clusterName BED
			my $BEDBlkName = $clusterName."|j=".$totalJNum."|read=".$prominentJunctSupportReadNum;
			
			my $BEDLine = junctInfoToBEDLine($prominentJunctStr, $NGSJunctStrndHsh{$prominentJunctStr}, 20, 20, $BEDBlkName, 10000);
			print JCLUSBED $BEDLine."\n";
			
		}
		close JCLUSBED;
		
		#---defining overlapping junction clusters
		print "Defining strand specific overlapping of clusters\n";
		my %jClusterSSOverlapHsh;
		foreach my $refJunctStr (sort {$a cmp $b} keys %SSOvrlpJunctIntronHsh) {
			my $refCluster = $jClusterNameByJunctHsh{$refJunctStr};
			foreach my $qryJunctStr (sort {$a cmp $b} keys %{$SSOvrlpJunctIntronHsh{$refJunctStr}}) {
				my $qryCluster = $jClusterNameByJunctHsh{$qryJunctStr};
				${$jClusterSSOverlapHsh{$refCluster}}{$qryCluster}++ if (($refCluster ne $qryCluster) and (${$jClusterInfoHsh{$qryCluster}}{"prominentJunct"} eq $qryJunctStr) and (${$jClusterInfoHsh{$refCluster}}{"prominentJunct"} eq $refJunctStr));
			}
		}

		my (%jClusterAntisenseOverlapHsh, %junctStrAntisenseOverlapHsh);
		
		print "Defining overlapping sense/Antisense jucntions clusters\n";
		
		foreach my $refJunctStr (keys %XSOvrlpJunctHsh) {
			my $refStrand = $NGSJunctStrndHsh{$refJunctStr};
			my $refCluster = $jClusterNameByJunctHsh{$refJunctStr};
			
			#---skip this refJunct if it is not strand specific
			next if (($refStrand ne "+") and  ($refStrand ne "-"));
			foreach my $qryJunctStr (keys %{$XSOvrlpJunctHsh{$refJunctStr}}) {
				my $qryStrand = $NGSJunctStrndHsh{$qryJunctStr};
				my $qryCluster = $jClusterNameByJunctHsh{$qryJunctStr};

				#---skip this refJunct if it is not strand specific
				next if (($qryStrand ne "+") and  ($qryStrand ne "-"));

				if (($refStrand ne $qryStrand) and ($refCluster ne $qryCluster)) {
					${$jClusterAntisenseOverlapHsh{$refCluster}}{$qryCluster}++;
					${$junctStrAntisenseOverlapHsh{$refJunctStr}}{$qryJunctStr}++;
				}
			}
		}

		#---print the overlapping antisense junctions
		open (TMPLOG, ">$outDir/ovrlpAntisenseJunct/overlappingAntisenseJunctions.txt");
		foreach my $refJunctStr (sort {$a cmp $b} keys %junctStrAntisenseOverlapHsh) {
			my $refCluster = $jClusterNameByJunctHsh{$refJunctStr};
			my $refReadNum = $junctReadNumHsh{$refJunctStr};
			my $refIsMajor = "n";
			$refIsMajor = "y" if (exists $majorJunctAvgSplcingRatioHsh{$refJunctStr});
			foreach my $qryJunctStr (sort {$a cmp $b} keys %{$junctStrAntisenseOverlapHsh{$refJunctStr}}) {
				my $qryCluster = $jClusterNameByJunctHsh{$qryJunctStr};
				my $qryReadNum = $junctReadNumHsh{$qryJunctStr};
				my $qryIsMajor = "n";
				$qryIsMajor = "y" if (exists $majorJunctAvgSplcingRatioHsh{$qryJunctStr});
				my $readRatio = sprintf "%.06f", $refReadNum/$qryReadNum;
				print TMPLOG join "\t", ($refJunctStr, $qryJunctStr, $refCluster, $qryCluster, $refReadNum, $qryReadNum, $refIsMajor, $qryIsMajor, $readRatio."\n");
			}
		}
		close TMPLOG;

		open (readNumRATIO, ">$outDir/correlation/inferiorClusterReadNumRatio.txt");

		my %NGSJunctInfoHsh;
		
		foreach my $clusterName (sort {$a cmp $b} keys %clusterAllJunctHsh) {
			my $clusterPJunctRefHit = my $clusterAllJunctRefHit = "n";
			my @tmpOverlappingSSClusterAry = ();
			if (exists $jClusterSSOverlapHsh{$clusterName}) {
				foreach my $overlappingSSCluster (keys %{$jClusterSSOverlapHsh{$clusterName}}) {
					push @tmpOverlappingSSClusterAry, $overlappingSSCluster;
				}
			} else {
				push @tmpOverlappingSSClusterAry, "null";
			}

			my @tmpOverlappingAntisenseClusterAry = ();
			if (exists $jClusterAntisenseOverlapHsh{$clusterName}) {
				foreach my $overlappingAntisenseCluster (keys %{$jClusterAntisenseOverlapHsh{$clusterName}}) {
					push @tmpOverlappingAntisenseClusterAry, $overlappingAntisenseCluster;
				}
			} else {
				push @tmpOverlappingAntisenseClusterAry, "null";
			}

			my $clstPrmntReadNum = $junctReadNumHsh{${$jClusterInfoHsh{$clusterName}}{"prominentJunct"}};
			my $clstJunctNum = @{$clusterAllJunctHsh{$clusterName}};
			
			#----consider the junctions
			my $superCluster = "n";
			my $totalRdNumInCluster = 0;
			
			foreach my $junctStr (@{$clusterAllJunctHsh{$clusterName}}) {
				${$NGSJunctInfoHsh{$junctStr}}{"prmntReadNumRatio"} = sprintf "%0.8f", ($junctReadNumHsh{$junctStr}/$clstPrmntReadNum);
				${$NGSJunctInfoHsh{$junctStr}}{"clstJunctNum"} = $clstJunctNum;
				${$NGSJunctInfoHsh{$junctStr}}{"readNum"} = $junctReadNumHsh{$junctStr};
				my ($junctCntg, $junctStart, $junctStop) = split /:/, $junctStr;
				my $intronSize = $junctStop - $junctStart + 1;
				my $intronFrame = '';
				
				if ($intronSize/3 == int ($intronSize/3)) {
					$intronFrame = '3n';

				} elsif (($intronSize+1)/3 == int (($intronSize+1)/3)) {
					$intronFrame = '3n+1';
				
				} elsif (($intronSize+2)/3 == int (($intronSize+2)/3)) {
					$intronFrame = '3n+2';
				
				} else {
					die "Impossible intronFrame\n";
				}
				
				${$NGSJunctInfoHsh{$junctStr}}{"intronSize"} = $intronSize;
				${$NGSJunctInfoHsh{$junctStr}}{"intronFrame"} = $intronFrame;
				
				$totalRdNumInCluster += ${$NGSJunctInfoHsh{$junctStr}}{"readNum"};
				
				${$NGSJunctInfoHsh{$junctStr}}{"score"} = $junctScoreHsh{$junctStr};
				${$NGSJunctInfoHsh{$junctStr}}{"strnd"} = $NGSJunctStrndHsh{$junctStr};
				${$NGSJunctInfoHsh{$junctStr}}{"leftFlankSeq"} = ${$NGSJunctSeqHsh{$junctStr}}{"leftFlankSeq"};
				${$NGSJunctInfoHsh{$junctStr}}{"rightFlankSeq"} = ${$NGSJunctSeqHsh{$junctStr}}{"rightFlankSeq"};
				${$NGSJunctInfoHsh{$junctStr}}{"completeSeq"} = ${$NGSJunctSeqHsh{$junctStr}}{"completeSeq"};
				${$NGSJunctInfoHsh{$junctStr}}{"upSite"} = ${$NGSJunctSeqHsh{$junctStr}}{"upSite"};
				${$NGSJunctInfoHsh{$junctStr}}{"downSite"} = ${$NGSJunctSeqHsh{$junctStr}}{"downSite"};
				${$NGSJunctInfoHsh{$junctStr}}{"unq"} = $NGSJunctUnqHsh{$junctStr};
				${$NGSJunctInfoHsh{$junctStr}}{"BEDCoverRng"} = $NGSJunctBEDCoverRngHsh{$junctStr};
				${$NGSJunctInfoHsh{$junctStr}}{"cluster"} = $clusterName;
				${$NGSJunctInfoHsh{$junctStr}}{"majorJunct"} = "n"; 
				${$NGSJunctInfoHsh{$junctStr}}{"prmntInCluster"} = "n"; 
				${$NGSJunctInfoHsh{$junctStr}}{"onPrmntIsofm"} = "n"; #---if both $majorJunct and $prmntInCluster are y
				${$NGSJunctInfoHsh{$junctStr}}{"superJunct"} = "n";
				${$NGSJunctInfoHsh{$junctStr}}{"hitRef"} = "null";
				${$NGSJunctInfoHsh{$junctStr}}{"majorCluster"} = ${$jClusterInfoHsh{$clusterName}}{"majorCluster"};#---check if the cluster is a major cluster
				${$NGSJunctInfoHsh{$junctStr}}{"splicingRatio"} = $allJunctAvgSplcingRatioHsh{$junctStr};
				${$NGSJunctInfoHsh{$junctStr}}{"majorJunct"} = "y" if (exists $majorJunctAvgSplcingRatioHsh{$junctStr});
				${$NGSJunctInfoHsh{$junctStr}}{"prmntInCluster"} = "y" if (${$jClusterInfoHsh{$clusterName}}{"prominentJunct"} eq $junctStr); 
				${$NGSJunctInfoHsh{$junctStr}}{"onPrmntIsofm"} = "y" if ((${$NGSJunctInfoHsh{$junctStr}}{"majorJunct"} eq "y") and (${$NGSJunctInfoHsh{$junctStr}}{"prmntInCluster"} eq "y")); 

				${$NGSJunctInfoHsh{$junctStr}}{"loCmPrptn"} = ${$NGSJunctSeqHsh{$junctStr}}{"loCmPrptn"};
				${$NGSJunctInfoHsh{$junctStr}}{"leftExonSeq"} = ${$NGSJunctSeqHsh{$junctStr}}{"leftExonSeq"};
				${$NGSJunctInfoHsh{$junctStr}}{"rightExonSeq"} = ${$NGSJunctSeqHsh{$junctStr}}{"rightExonSeq"};
				${$NGSJunctInfoHsh{$junctStr}}{"leftExonPolyBasePrptn"} = ${$NGSJunctSeqHsh{$junctStr}}{"leftExonPolyBasePrptn"};
				${$NGSJunctInfoHsh{$junctStr}}{"rightExonPolyBasePrptn"} = ${$NGSJunctSeqHsh{$junctStr}}{"rightExonPolyBasePrptn"};

				if (exists $superJunctHsh{$junctStr}) {
					${$NGSJunctInfoHsh{$junctStr}}{"superJunct"} = "y";
					$superCluster = "y";
				}
				${$NGSJunctInfoHsh{$junctStr}}{"prmntDiff5"} = "null";
				${$NGSJunctInfoHsh{$junctStr}}{"prmntDiff3"} = "null";
				${$NGSJunctInfoHsh{$junctStr}}{"prmntDiffBoth"} = "null";

				if (exists $refGeneIDByJunctStrHsh{$junctStr}) {#----same junctStr in refGff
					${$NGSJunctInfoHsh{$junctStr}}{"hitRef"} = ${$jClusterInfoHsh{$clusterName}}{"hitRef"} = $refGeneIDByJunctStrHsh{$junctStr};
				}
				
				if ((${$NGSJunctInfoHsh{$junctStr}}{"prmntInCluster"} eq "n") and ($clstJunctNum > 1) and (${$NGSJunctInfoHsh{$junctStr}}{"superJunct"} eq "n") and (${$NGSJunctInfoHsh{$junctStr}}{"majorCluster"} eq "y")) {
					print readNumRATIO ${$NGSJunctInfoHsh{$junctStr}}{"prmntReadNumRatio"}."\n";
				}
			}

			@{${$jClusterInfoHsh{$clusterName}}{"allJunctStr"}} = @{$clusterAllJunctHsh{$clusterName}};
			@{${$jClusterInfoHsh{$clusterName}}{"ovrlpAS"}} = @tmpOverlappingAntisenseClusterAry;
			@{${$jClusterInfoHsh{$clusterName}}{"ovrlpSS"}} = @tmpOverlappingSSClusterAry;
			${$jClusterInfoHsh{$clusterName}}{"unq"} = $NGSJunctUnqHsh{${$jClusterInfoHsh{$clusterName}}{"prominentJunct"}};
			${$jClusterInfoHsh{$clusterName}}{"hitRef"} = "null";
			${$jClusterInfoHsh{$clusterName}}{"superCluster"} = $superCluster;
			${$jClusterInfoHsh{$clusterName}}{"totalRdNumInCluster"} = $totalRdNumInCluster;

		}
		close readNumRATIO;

		return \%clusterAllJunctHsh, \%jClusterNameByJunctHsh, \%jClusterInfoHsh, \%jClusterSSOverlapHsh, \%NGSJunctInfoHsh;
	}
	############################################## findSuperJunctionOverlapCluster ########################################################
	sub defineExonSkipping {
	
		print "Defining exon skipping events.\n";
	
		my %superJunctHsh = %{$_[0]};
		my %majorJunctAvgSplcingRatioHsh = %{$_[1]};
		my %jClusterNameByJunctHsh = %{$_[2]};
		
		#---hash to return
		my (%totalExonSkippingTypeHsh, %majorExonSkippingTypeHsh, %exactExonSkippingClusterHsh, %superJunctOvrlpClusterHsh);
		
		foreach my $superJunct (sort {$a cmp $b} keys %superJunctHsh) {
			#--- non-super junctions are defined as inferior junctions here
			my $startTotalHit = my $endTotalHit = my $startMajorHit = my $endMajorHit = 0;
			my @superJunctSplt = split /:/, $superJunct;
			my $superJunctStart = $superJunctSplt[1];
			my $superJunctEnd = $superJunctSplt[2];;

			#---check for exact overlapping of start and end of the inferior junctions
			foreach my $inferiorJunctStr (sort {$a cmp $b} keys %{$SSOvrlpJunctIntronHsh{$superJunct}}) {
				my @inferiorJunctStrSplt = split /:/, $inferiorJunctStr;
				my $inferiorJunctStart = $inferiorJunctStrSplt[1];
				my $inferiorJunctEnd = $inferiorJunctStrSplt[2];
					
				my $inferiorClusterNum = $jClusterNameByJunctHsh{$inferiorJunctStr};
					 
				${$superJunctOvrlpClusterHsh{$superJunct}}{$inferiorClusterNum}++;
					
				if ($superJunctStart eq $inferiorJunctStart) {
					${$exactExonSkippingClusterHsh{$superJunct}}{$inferiorJunctStr}++;
					$startTotalHit ++;
					$startMajorHit ++ if (exists $majorJunctAvgSplcingRatioHsh{$inferiorJunctStr});
				}
				
				if ($superJunctEnd eq $inferiorJunctEnd) {
					${$exactExonSkippingClusterHsh{$superJunct}}{$inferiorJunctStr}++;
					$endTotalHit ++;
					$endMajorHit ++ if (exists $majorJunctAvgSplcingRatioHsh{$inferiorJunctStr});
				}
			} #---end of foreach my $inferiorJunctStr (sort {$a cmp $b} keys %{$SSOvrlpJunctIntronHsh{$superJunct}}) {
			
			if (($startTotalHit > 0) and ($endTotalHit > 0)) {#---exact exonSkipping
				$totalExonSkippingTypeHsh{$superJunct} = "exact";

			} elsif (($startTotalHit > 0) or ($endTotalHit > 0)) {#---partial exonSkpping
				$totalExonSkippingTypeHsh{$superJunct} = "partial";
			
			} elsif (($startTotalHit == 0) or ($endTotalHit  == 0)) {#---inexact exonSkpping
				$totalExonSkippingTypeHsh{$superJunct} = "inexact";
				
			} else {#---unexpected events
				die "Unexpected sincerio of exon skipping events.\n";
			}
			
			if (($startMajorHit > 0) and ($endMajorHit > 0)) {#---exact exonSkipping
				$majorExonSkippingTypeHsh{$superJunct} = "exact";

			} elsif (($startMajorHit > 0) or ($endMajorHit > 0)) {#---partial exonSkpping
				$majorExonSkippingTypeHsh{$superJunct} = "partial";
			
			} elsif (($startMajorHit == 0) or ($endMajorHit  == 0)) {#---inexact exonSkpping
				$majorExonSkippingTypeHsh{$superJunct} = "inexact";
				
			} else {#---unexpected events
				die "Unexpected sincerio of exon skipping events.\n";
			}

		}
		
		open (TMPLOG, ">$outDir/exonSkip/totalExonSkippingEvents.txt");
		foreach my $superJunct (sort {$totalExonSkippingTypeHsh{$a} cmp $totalExonSkippingTypeHsh{$b}} keys %totalExonSkippingTypeHsh) {
			print TMPLOG $superJunct."\t".$totalExonSkippingTypeHsh{$superJunct};
			my $skippedClusterNum = keys %{$superJunctOvrlpClusterHsh{$superJunct}};
			foreach my $inferiorClusterNum (sort {$a cmp $b} keys %{$superJunctOvrlpClusterHsh{$superJunct}}) {
				print TMPLOG "\t".$inferiorClusterNum;
			}
			print TMPLOG "\n";
		}
		close TMPLOG;

		open (TMPLOG, ">$outDir/exonSkip/majorExonSkippingEvents.txt");
		foreach my $superJunct (sort {$majorExonSkippingTypeHsh{$a} cmp $majorExonSkippingTypeHsh{$b}} keys %majorExonSkippingTypeHsh) {
			print TMPLOG $superJunct."\t".$majorExonSkippingTypeHsh{$superJunct};
			my $skippedClusterNum = keys %{$superJunctOvrlpClusterHsh{$superJunct}};
			foreach my $inferiorClusterNum (sort {$a cmp $b} keys %{$superJunctOvrlpClusterHsh{$superJunct}}) {
				print TMPLOG "\t".$inferiorClusterNum;
			}
			print TMPLOG "\n";
		}
		close TMPLOG;

		return (\%totalExonSkippingTypeHsh, \%majorExonSkippingTypeHsh, \%exactExonSkippingClusterHsh, \%superJunctOvrlpClusterHsh);
		
	}
	############################################## defineAltSplicingSiteWithinCluster ########################################################
	sub defineAltSplicingSiteWithinCluster {
	
		my %clusterAllJunctHsh = %{$_[0]};
		my %majorJunctAvgSplcingRatioHsh = %{$_[1]};
		my %NGSJunctInfoHsh = %{$_[2]};
		
		print "Defining alternative splicing site within intron clusters\n";
		
		#---hashes to return
		my (%splicingSiteDiffAllHsh, %splicingSiteDiffHshByClusterHsh, %splicingSiteDiffMajorHsh);

		foreach my $clusterName (sort {$a cmp $b} keys %clusterAllJunctHsh) {
			
			if (@{$clusterAllJunctHsh{$clusterName}} > 1) { #---clusters with more than 1 junct, not orphan and superJunctions
				my $lastIndex = $#{$clusterAllJunctHsh{$clusterName}};
				
				#---store the Strnd of the first junction, to make sure all junct within cluster is on the same strd
				#---since 0 wont appear as $i,so get the index 0 strd first
				my $clusterStrnd = $NGSJunctStrndHsh{${$clusterAllJunctHsh{$clusterName}}[0]};
				
				#---pairwise comparison of cluster junctions
				for (my $i=$lastIndex; $i>0; $i--) {#---5, 4, 3, 2, 1 for 6 junctions
					my $refJunctStr = ${$clusterAllJunctHsh{$clusterName}}[$i];
					my @refJunctStrSplt = split /:/, $refJunctStr;
					my $refJunctStart = $refJunctStrSplt[1];
					my $refJunctEnd = $refJunctStrSplt[2];
					
					my $refIsMajor = "n";
					$refIsMajor = "y" if (exists $majorJunctAvgSplcingRatioHsh{$refJunctStr});
					
					my $refIsPrmntInClust = ${$NGSJunctInfoHsh{$refJunctStr}}{"prmntInCluster"};
					
					die "Intron cluster contain jucntions of diifferent strands. Unexpected scinerio.\n" if ($clusterStrnd ne $NGSJunctStrndHsh{$refJunctStr});

					for (my $j=$i-1; $j>=0; $j--) {#----4, 3, 2, 1, 0|| 3, 2, 1, 0 || 2, 1, 0 || 1, 0|| 0 for 6 junctions
						my $qryJunctStr = ${$clusterAllJunctHsh{$clusterName}}[$j];
						my @qryJunctStrSplt = split /:/, $qryJunctStr;
						my $qryJunctStart = $qryJunctStrSplt[1];
						my $qryJunctEnd = $qryJunctStrSplt[2];
						
						my $startDiff = ($refJunctStart - $qryJunctStart);
						my $endDiff = ($qryJunctEnd - $refJunctEnd);
						my $fullDiff = $endDiff + $startDiff;
						${$splicingSiteDiffAllHsh{"full"}}{$fullDiff}++;
						${${$splicingSiteDiffHshByClusterHsh{$clusterName}}{"full"}}{$fullDiff}++;

						if ($refIsMajor eq "y") {
							${$splicingSiteDiffMajorHsh{"full"}}{$fullDiff}++;
						}

						if ($clusterStrnd eq "-") {#---for minus strand
		
							${$splicingSiteDiffAllHsh{"3"}}{$startDiff}++;
							${$splicingSiteDiffAllHsh{"5"}}{$endDiff}++;
							${${$splicingSiteDiffHshByClusterHsh{$clusterName}}{"3"}}{$startDiff}++;
							${${$splicingSiteDiffHshByClusterHsh{$clusterName}}{"5"}}{$endDiff}++;
							
							if ($refIsMajor eq "y") {
								${$splicingSiteDiffMajorHsh{"3"}}{$startDiff}++;
								${$splicingSiteDiffMajorHsh{"5"}}{$endDiff}++;							
							}
							
							if ($refIsPrmntInClust eq "y") {
								${$NGSJunctInfoHsh{$qryJunctStr}}{"prmntDiff5"} = $endDiff;
								${$NGSJunctInfoHsh{$qryJunctStr}}{"prmntDiff3"} = $startDiff;
								${$NGSJunctInfoHsh{$qryJunctStr}}{"prmntDiffBoth"} = $fullDiff;
							}
							
						} else {#---for both no strand and plus strand

							${$splicingSiteDiffAllHsh{"5"}}{$startDiff}++;
							${$splicingSiteDiffAllHsh{"3"}}{$endDiff}++;
							${${$splicingSiteDiffHshByClusterHsh{$clusterName}}{"5"}}{$startDiff}++;
							${${$splicingSiteDiffHshByClusterHsh{$clusterName}}{"3"}}{$endDiff}++;
							
							if ($refIsMajor eq "y") {
								${$splicingSiteDiffMajorHsh{"5"}}{$startDiff}++;
								${$splicingSiteDiffMajorHsh{"3"}}{$endDiff}++;					
							}
							
							if ($refIsPrmntInClust eq "y") {
								${$NGSJunctInfoHsh{$qryJunctStr}}{"prmntDiff5"} = $startDiff;
								${$NGSJunctInfoHsh{$qryJunctStr}}{"prmntDiff3"} = $endDiff;
								${$NGSJunctInfoHsh{$qryJunctStr}}{"prmntDiffBoth"} = $fullDiff;
							}

						}
					}
				}
			}#---end of if (@{$clusterAllJunctHsh{$clusterName}} > 1)
		}#---end of foreach my $clusterName (sort {$a cmp $b} keys %clusterAllJunctHsh) {
=pod
		#---plot the bar chart for all
		my $plotCutoff = 30;
		foreach my $bound53 (keys %splicingSiteDiffAllHsh) {
			
			my (%tmpForPlotFullHsh, %tmpForPlotCutoffHsh);
			
			foreach my $diff (keys  %{$splicingSiteDiffAllHsh{$bound53}}) {
				
				$tmpForPlotFullHsh{$diff} = ${$splicingSiteDiffAllHsh{$bound53}}{$diff};
				
				if (($diff <= $plotCutoff) and ($diff > 0)) {
					$tmpForPlotCutoffHsh{$diff} = ${$splicingSiteDiffAllHsh{$bound53}}{$diff};
				} elsif ($diff > $plotCutoff) {
					$tmpForPlotCutoffHsh{$plotCutoff} = 0 if (not exists $tmpForPlotCutoffHsh{$plotCutoff});
					$tmpForPlotCutoffHsh{$plotCutoff} += ${$splicingSiteDiffAllHsh{$bound53}}{$diff};
				}
			}
			
			GNUPlotBarChartNumberItem(\%tmpForPlotFullHsh, "intron $bound53 end altSplicingSiteDistance", "$outDir/plotData/intron$bound53.altSplicingSiteDistance.full.dat", "$outDir/plotPdf/intron$bound53.altSplicingSiteDistance.full.pdf", "Distance (nt)", "Frequency", "n", "n" );
			GNUPlotBarChartNumberItem(\%tmpForPlotCutoffHsh, "intron $bound53 end altSplicingSiteDistance", "$outDir/plotData/intron$bound53.altSplicingSiteDistance.cutoff.dat", "$outDir/plotPdf/intron$bound53.altSplicingSiteDistance.cutoff.pdf", "Distance (nt)", "Frequency", "n", "n");
		}

		#---plot the bar chart for major
		foreach my $bound53 (keys %splicingSiteDiffMajorHsh) {
			
			my (%tmpForPlotFullHsh, %tmpForPlotCutoffHsh);
			
			foreach my $diff (keys  %{$splicingSiteDiffMajorHsh{$bound53}}) {
				
				$tmpForPlotFullHsh{$diff} = ${$splicingSiteDiffMajorHsh{$bound53}}{$diff};
				
				if (($diff <= $plotCutoff) and ($diff > 0)) {
					$tmpForPlotCutoffHsh{$diff} = ${$splicingSiteDiffMajorHsh{$bound53}}{$diff};
				} elsif ($diff > $plotCutoff) {
					$tmpForPlotCutoffHsh{$plotCutoff} = 0 if (not exists $tmpForPlotCutoffHsh{$plotCutoff});
					$tmpForPlotCutoffHsh{$plotCutoff} += ${$splicingSiteDiffMajorHsh{$bound53}}{$diff};
				}
			}
			
			GNUPlotBarChartNumberItem(\%tmpForPlotFullHsh, "major intron $bound53 end altSplicingSiteDistance", "$outDir/plotData/intron$bound53.major.altSplicingSiteDistance.full.dat", "$outDir/plotPdf/intron$bound53.major.altSplicingSiteDistance.full.pdf", "Distance (nt)", "Frequency", "n", "n");
			GNUPlotBarChartNumberItem(\%tmpForPlotCutoffHsh, "major intron $bound53 end altSplicingSiteDistance", "$outDir/plotData/intron$bound53.major.altSplicingSiteDistance.cutoff.dat", "$outDir/plotPdf/intron$bound53.major.altSplicingSiteDistance.cutoff.pdf", "Distance (nt)", "Frequency", "n", "n");
		}
=cut
		return (\%splicingSiteDiffAllHsh, \%splicingSiteDiffHshByClusterHsh, \%NGSJunctInfoHsh);
	
	}
	############################################## defineJunctionsOfProminentIsoform ########################################################
	sub defineJunctionsOfProminentIsoform {#---junctions on a prominent isoform is defined as a junction that has the highest number of read among all overlapping junctions

		#---rationale:
		#---check for each superjunction, to check whether it has the highest read num among all overlapping clusters.
		#---if yes, the so_ and sc_cluster (i.e. superjunction) will be stored in the prominent cluster hash
		#---if no, the overlapping clusters will be stored in the prominent cluster hash
		#---then go through each cluster and pick up the junctions with the highest number of reads (i.e. prominent read)
		#---go through each prominent read to see if they overlap with each other, it may happen for two s_ clusters since superjunction clusters 
		
		print "Defining the junctions of the most prominent isoform.\n";
		my (%junctOnProminentIsofmHsh);
		
		my %jClusterNameByJunctHsh = %{$_[0]};
		my %jClusterInfoHsh = %{$_[1]};
		my %jClusterSSOverlapHsh = %{$_[2]};
		
		#----scan the superJunctionFirst;
		foreach my $clusterName (keys %jClusterInfoHsh) {
			if ($clusterName =~ m/^s/) {#---superJunction
				my $superClusterName = $clusterName;
				my $superSupportReadNum = ${$jClusterInfoHsh{$superClusterName}}{"mostSupportReadNum"};
				my $superIsProminent = "y";
				
				#----check the inferior clusters
				foreach my $overlapInferiorCluster (keys %{$jClusterSSOverlapHsh{$superClusterName}}) {
					my $inferiorSupportReadNum = ${$jClusterInfoHsh{$overlapInferiorCluster}}{"mostSupportReadNum"};
					if ($inferiorSupportReadNum > $superSupportReadNum) {
						$superIsProminent = "n";
					}
				}
				
				if ($superIsProminent eq "y") {
					#---skip the inferiorClsuter if the superJunction cluster is prominent;
					foreach my $overlapInferiorCluster (keys %{$jClusterSSOverlapHsh{$superClusterName}}) {
						${$jClusterInfoHsh{$overlapInferiorCluster}}{"onProminentIsoform"} = "n";#---to skip the inferior junct later on since the super junction is major
					}
					
					my $clusterIsMajor = ${$jClusterInfoHsh{$superClusterName}}{"majorCluster"};
					
					if ($clusterIsMajor eq "y") {
						my $prominentJunctStr = ${$jClusterInfoHsh{$superClusterName}}{"prominentJunct"};
						${$jClusterInfoHsh{$superClusterName}}{"onProminentIsoform"} = "y";
						$junctOnProminentIsofmHsh{$prominentJunctStr} = $superSupportReadNum;
						
					} else {#---superJunction has the major support num but not major
						${$jClusterInfoHsh{$superClusterName}}{"onProminentIsoform"} = "n";
					}
					
				} else { #---superJunction is not prominent
					${$jClusterInfoHsh{$superClusterName}}{"onProminentIsoform"} = "n";
				}
			}
		} #---end of foreach my $clusterName (keys %jClusterInfoHsh)
		
		#----scan the unprocessed inferior junctions i.e. not exists ${$jClusterInfoHsh{$clusterName}}{"onProminentIsoform"}
		foreach my $clusterName (keys %jClusterInfoHsh) {
			if (not exists ${$jClusterInfoHsh{$clusterName}}{"onProminentIsoform"}) {
				my $clusterIsMajor = ${$jClusterInfoHsh{$clusterName}}{"majorCluster"};
				if ($clusterIsMajor eq "y") {
					my $prominentJunctStr = ${$jClusterInfoHsh{$clusterName}}{"prominentJunct"};
					${$jClusterInfoHsh{$clusterName}}{"onProminentIsoform"} = "y";
					my $supportReadNum = ${$jClusterInfoHsh{$clusterName}}{"mostSupportReadNum"};
					$junctOnProminentIsofmHsh{$prominentJunctStr} = $supportReadNum;
				} else {
					${$jClusterInfoHsh{$clusterName}}{"onProminentIsoform"} = "n";
				}
			}
		}
		
		#---Just double check to make sure the junctions on the prominent isofrom doesnt overlapp with each other
		print "Making sure the junctions on the prominent isofrom do not overlap with each other.\n";
		my %tmpCheckJunctHsh;
		foreach my $refJunctStr (keys %junctOnProminentIsofmHsh) {
			$tmpCheckJunctHsh{$refJunctStr}++;
			foreach my $qryJunctStr (keys %junctOnProminentIsofmHsh) {
				next if (exists $tmpCheckJunctHsh{$qryJunctStr});
				if ($refJunctStr ne $qryJunctStr) {
					die "Unexpected scenerio: Junctions on prominent isoform overlap with each other.\n" if (exists ${$SSOvrlpJunctIntronHsh{$refJunctStr}}{$qryJunctStr});
				}
			}
		}
		
		%tmpCheckJunctHsh = ();
		
		#---print the junction on prominent isoform
		open (TMPLOG, ">$outDir/junctInfo/junctOnProminentIsofm.txt");
		foreach my $junctStr (sort{$jClusterNameByJunctHsh{$a} cmp $jClusterNameByJunctHsh{$b}} keys %junctOnProminentIsofmHsh) {
			print TMPLOG $junctStr."\t".$jClusterNameByJunctHsh{$junctStr}."\n";
		}
		close TMPLOG;
		
		return (\%jClusterInfoHsh, \%junctOnProminentIsofmHsh);
		
	}
	
}
########################################################################## GNUPlotXYScatterWithLines
sub GNUPlotXYScatterWithLines {

	my %XYHsh = %{$_[0]};
	my $plotFilePath = $_[1];
	my $plotDataPath = $_[2];
	my $xlable = $_[3];
	my $ylable = $_[4];
	my $xscale = $_[5];
	my $yscale = $_[6];
	my $title = $_[7];
	my $verbose = $_[8];
	
	my $GNULogXCmd = "";
	$GNULogXCmd = "set logscale x" if ($xscale eq "log");
	my $GNULogYCmd = "";
	$GNULogYCmd = "set logscale y" if ($yscale eq "log");

	$plotFilePath .= ".pdf" if ($plotFilePath !~ m/\.pdf$/);

	my @filePathSplt = split /\//, $plotFilePath;
	my $fileName = $filePathSplt[-1];

	print "Running GNUPlotXYScatterWithLines for $fileName.\n" if ($verbose eq "y");
	
	#---creat a tmp file
	open (TMPFILE, ">$plotDataPath");
	for my $x (sort {$a <=> $b} keys %XYHsh) {
		print TMPFILE $x."\t".$XYHsh{$x}."\n";
	}
	close TMPFILE;
	
	#---do the GNUPLOT
	open (GNUPLOT, "|gnuplot");
	print GNUPLOT <<EOPLOT;
	set terminal postscript color solid
	set output "| ps2pdf - $plotFilePath 2>/dev/null";
	unset logscale x; 
	unset logscale y; 
	$GNULogXCmd;
	$GNULogYCmd;
	set xlabel "$xlable";
	set ylabel "$ylable";
	set title "$title";
	set nokey;
   	plot '$plotDataPath' using 1:2 with lines;
EOPLOT
	close(GNUPLOT);
	#rmtree(['tmp.dat'], 0, 1); #---non-verbose removal of tmp file
}
########################################################################## GNUPLOTAryHistogram
sub GNUPLOTAryHistogram {

	my @numAry = @{$_[0]};
	my $plotFilePath = $_[1];
	my $plotDataPath = $_[2];
	my $xlable = $_[3];
	my $ylable = $_[4];
	my $title = $_[5];
	my $verbose = $_[6];
	
	my @filePathSplt = split /\//, $plotFilePath;
	my $fileName = $filePathSplt[-1];
	
	print "Running GNUPLOTAryHistogram for $fileName.\n" if ($verbose eq "y");
	
	#---calculate the optimal bin
	@numAry = sort {$a <=> $b} @numAry;
	my $tail5PctPos = int ($#numAry*0.05);
	my $upper95PctPos = $#numAry - $tail5PctPos;
	my $lower95PctPos = $tail5PctPos;
	my $binWidth = ($numAry[$upper95PctPos] - $numAry[$lower95PctPos])/100;
	#$binWidth = 1 if ($binWidth < 1);

	#---creat a tmp file
	open (TMPFILE, ">$plotDataPath");
	for (@numAry) {print TMPFILE $_."\n";}
	close TMPFILE;
	
	#---do the GNUPLOT
	open (GNUPLOT, "|gnuplot");
	print GNUPLOT <<EOPLOT;
	set terminal postscript color solid
	bin_width = $binWidth
	bin_number(x) = floor(x/bin_width)
	rounded(x) = bin_width * (bin_number(x) + 0.5)
	UNITY = 1
	unset logscale x;
	unset logscale y;
	set output "| ps2pdf - $plotFilePath 2>/dev/null";
	set xlabel "$xlable";
	set ylabel "$ylable";
	set title "$title";
   	plot '$plotDataPath' u (rounded(\$1)):(UNITY) t 'bin at $binWidth' smooth frequency w histeps
EOPLOT
	close(GNUPLOT);
	#rmtree(['./tmpRand.dat'], 0, 1); #---non-verbose removal of tmp file
}
########################################################################## generateCumulativePctHash
sub generateCumulativePctHash {

	#---both the keys and the values have to numbers
	my %countHsh = %{$_[0]};
	my $order = $_[1]; #---a for ascending or d for descending

	#---find the total
	my $totalCount = 0;
	foreach my $item (sort {$a <=> $b} keys %countHsh) {$totalCount += $countHsh{$item};}
	
	#---generate the cumulative values;
	my $cumulativeCount = 0;
	my %cumulativeHsh;
	if ($order eq "a") {
		foreach my $item (sort {$a <=> $b} keys %countHsh) {
			$cumulativeCount += $countHsh{$item};
			my $cumulativePct = sprintf "%.06f", ($cumulativeCount/$totalCount)*100;
			$cumulativeHsh{$item} = $cumulativePct;
		}

	} elsif ($order eq "d") {
		foreach my $item (sort {$b <=> $a} keys %countHsh) {
			$cumulativeCount += $countHsh{$item};
			my $cumulativePct = sprintf "%.06f", ($cumulativeCount/$totalCount)*100;
			$cumulativeHsh{$item} = $cumulativePct;
		}
	}	
	
	return \%cumulativeHsh;
}
########################################################################## GNUPlotBarChartNumberItem
sub GNUPlotBarChartNumberItem {

	#---assuming the keys of the hash are numbers
	my %barDataHsh = %{$_[0]};
	my $plotTitle = $_[1];
	my $dataFilePath = $_[2];
	my $plotFilePath = $_[3];
	my $xLabel = $_[4];
	my $yLabel = $_[5];
	my $removeDataFile = $_[6];
	my $verbose = $_[7];
	
	#---add extension to =plot filename if it doesnt have
	$plotFilePath .= ".pdf" if ($plotFilePath !~ m/\.pdf$/);
	
	my @filePathSplt = split /\//, $plotFilePath;
	my $fileName = $filePathSplt[-1];

	print "Running GNUPlotBarChartNumberItem for $fileName.\n" if ($verbose eq "y");

	#---get the total count
	my $totalCount = 0;
	foreach my $item (sort {$a <=> $b} keys %barDataHsh) {$totalCount += $barDataHsh{$item};}

	#---generate a tmpData file
	my $tmpBarDatPath = $dataFilePath;
	open (TMPBARDATA, ">$tmpBarDatPath");
	my $maxPct = 0;
	foreach my $item (sort {$a <=> $b} keys %barDataHsh) {
		my $itemCount = $barDataHsh{$item};
		my $tmpPct = sprintf "%.2f", ($itemCount/$totalCount)*100;
		print TMPBARDATA $item."\t".$tmpPct."\t".$itemCount."\n";
		$maxPct = $tmpPct if ($tmpPct > $maxPct);
	}
	close TMPBARDATA;
	
	#---define the max y values
	my $yMax = ((int ($maxPct/10))*10) + 10; #---in scale of 10%
	my $ytics = int ($yMax/10);
	my $y2Max = int ($totalCount*($yMax/100));
	my $y2tics = int ($y2Max/10);

	open (GNUPLOT, "|gnuplot");
	print GNUPLOT <<EOPLOT;
	set terminal postscript color solid
	set output "| ps2pdf - $plotFilePath 2>/dev/null";
	set key left;
	set title "$plotTitle";
	set xlabel "$xLabel";
	set ylabel "% $yLabel";
	set y2label "count $yLabel";
	set yrange [0:$yMax];
	set ytics $ytics nomirror;
	set y2range [0:$y2Max];
	set y2tics $y2tics nomirror;
	set boxwidth 1 absolute;
	set style fill solid 2.00 border 0;
	set style histogram;
	set style data histograms;
	plot '$tmpBarDatPath' using 2:xtic(1);
EOPLOT
	
	system ("rm $tmpBarDatPath") if $removeDataFile eq "y";
	
}
########################################################################## getIntronFromExonRng
sub getIntronFromExonRng {

	#---incoming variables
	my %exonRngByGeneHsh = %{$_[0]};
	my %cntgByGeneHsh = %{$_[1]};
	my %strndByGeneHsh = %{$_[2]};
	my $strToPrint = $_[3];
	
	#---outgoing variables
	my (%intronRngByGeneHsh, %junctStrByIntronIDHsh, %XStrndRngByCntgByIntronID, %strndByIntronIDHsh, %geneIDByJunctStrHsh);
	
	#---on screen progress scale
	printProgressScale("\n$strToPrint", 50, 10);
	my $totalGeneNum = keys %exonRngByGeneHsh;
	my $progressGeneNum = 0;
	
	#---go through each gene, see if there's intron, and store the ranges
	foreach my $geneID (sort {$a cmp $b} keys %exonRngByGeneHsh) {

		#---update the on-screen progress
		$progressGeneNum++;
		updateProgressBar("$geneID", $progressGeneNum, $totalGeneNum, 50, 10);
		
		#---get the cntg and strnd
		my $cntg = $cntgByGeneHsh{$geneID};
		my $strnd = $strndByGeneHsh{$geneID};
		
		#---get all exon bounds
		my @tmpBoundAry;
		foreach my $exon (keys %{$exonRngByGeneHsh{$geneID}}) {
			push @tmpBoundAry, ${${$exonRngByGeneHsh{$geneID}}{$exon}}{"start"};
			push @tmpBoundAry, ${${$exonRngByGeneHsh{$geneID}}{$exon}}{"end"};
		}
		
		#---if more than one exon
		if (@tmpBoundAry > 2) { 
			#---sort the bounds
			my @sortedTmpBoundAry = sort {$a <=> $b} @tmpBoundAry;
			#---get the intronBounds, go to all odd number indexes and take itself and the +1 index values, i.e. $i = 1, 3, 5 if there're 0,1,2,3,4,5,6,7  
			my $intronNum = 0;
			for (my $i=1; $i<(@sortedTmpBoundAry-1); $i=$i+2) {
				$intronNum++;
				my $intronStart = $sortedTmpBoundAry[$i]+1;
				my $intronEnd = $sortedTmpBoundAry[$i+1]-1;
				my $junctStr = join ":", ($cntg, $intronStart, $intronEnd);
				#my $intronID = $geneID.":".$intronNum;
				my $intronID = $junctStr; #----ad hoc switch intron ID to junctStr
				${${$intronRngByGeneHsh{$geneID}}{$intronID}}{"start"} = $intronStart;
				${${$intronRngByGeneHsh{$geneID}}{$intronID}}{"end"} = $intronEnd;
				${${$XStrndRngByCntgByIntronID{$cntg}}{$intronID}}{"start"} = $intronStart;
				${${$XStrndRngByCntgByIntronID{$cntg}}{$intronID}}{"end"} = $intronEnd;
				$strndByIntronIDHsh{$intronID} = $strnd;
				$junctStrByIntronIDHsh{$intronID} = $junctStr;
				$geneIDByJunctStrHsh{$junctStr} = $geneID;
			}
		}
	}

	print "\n\n";

	return (\%intronRngByGeneHsh, \%junctStrByIntronIDHsh, \%XStrndRngByCntgByIntronID, \%strndByIntronIDHsh, \%geneIDByJunctStrHsh)
	
}
########################################################################## updateProgressBar
sub updateProgressBar {
	
	my $strToPrint = $_[0];
	my $progressCount = $_[1];
	my $totalCount = $_[2];
	my $scaleMax = $_[3];
	my $extraWhiteSpaceNum = $_[4]; #---erase the longer infos during the progress
	
	my $progressPct = int (($progressCount/$totalCount)*$scaleMax);

	my $progressBar = "|";
	for my $i (1..$progressPct) {$progressBar .= ">";}
	for my $i (($progressPct+1)..$scaleMax) {$progressBar .= " ";}
	$progressBar .= "|";

	my $extraWhiteSpaceStr = "";
	for my $i (1..$extraWhiteSpaceNum) {$extraWhiteSpaceStr .= " ";}
	
	print $progressBar.$strToPrint.$extraWhiteSpaceStr."\r";

}
########################################################################## printProgressScale
sub printProgressScale {

	my $strToPrint = $_[0];
	my $scaleMax = $_[1];

	my $scaleSpace = "|";
	for my $i (1..$scaleMax) {$scaleSpace .= "-";}
	$scaleSpace .= "|100%";
	
	print $strToPrint."\n";
	print $scaleSpace."\n";
}
########################################################################## summarizeAlternativeSplicingOnRefTranscript
sub summarizeAlternativeSplicingOnRefTranscript {

	my %refIntronRngXS = %{$_[0]}; #---${${$XStrndRngByCntgByIntronID{$cntg}}{$intronID}}{"start"} = $intronStart;
	my %refRngXSHsh = %{$_[1]}; #---${${$XStrndRngByCntgByWholeGeneHsh{$seq}}{$geneID}}{"start"} = $featureStart
	my %refStrndHsh = %{$_[2]}; #---$strndByGeneHsh{$geneID} = $geneStrnd;
	my %junctScoreHsh = %{$_[3]}; #---$junctScoreHsh{$junctStr} = $score;
	my %junctReadNumHsh = %{$_[4]}; #---$junctReadNumHsh{$junctStr} = $readNum;
	my %NGSJunctCntgHsh = %{$_[5]}; #---$junctCtngHsh{$junctStr} = $cntg;
	my %SSOvrlpRefJNGSJHsh = %{$_[6]}; #---${$SSHitByRefHsh{$refFtur}}{$qryFtur} = 0;
	my %XSOvrlpRefJNGSJHsh = %{$_[7]}; #---${$XStrndHitByRefHsh{$refFtur}}{$qryFtur} = 0;
	my %SSOvrlpRefTrnscptNGSJHsh = %{$_[8]}; #---${$SSHitByRefHsh{$refFtur}}{$qryFtur} = 0;
	my %XSOvrlpRefTrnscptNGSJHsh = %{$_[9]}; #---${$XStrndHitByRefHsh{$refFtur}}{$qryFtur} = 0;
	my %superJunctHsh = %{$_[10]}; #---$superJunctHsh{$refJunctStr}++ if ((not exists ${$SSOvrlpJunctIntronHsh{$refOverlapJunct}}{$qryJunctStr}) and ($qryJunctStr ne $refOverlapJunct));
	my %clusterAllJunctHsh = %{$_[11]}; #---push @{$clusterAllJunctHsh{$clusterName}}, $junctStr;
	my %jClusterNameByJunctHsh = %{$_[12]}; #---$jClusterNameByJunctHsh{$junctStr}  = $clusterName;
	my %jClusterInfoHsh = %{$_[13]}; #---${$jClusterInfoHsh{$clusterName}}{"prominentJunct"} = $prominentJunctStr; etc
	my %totalExonSkippingTypeHsh = %{$_[14]}; #---$totalExonSkippingTypeHsh{$superJunct} = "exact";
	my %superJunctOvrlpClusterHsh = %{$_[15]}; #---${$superJunctOvrlpClusterHsh{$superJunct}}{$inferiorClusterNum}++;
	my %SSOvrlpJunctIntronHsh = %{$_[16]}; #---${$SSHitByRefHsh{$refFtur}}{$qryFtur} = 0;
	my %refExonRngByGeneHsh = %{$_[17]}; #---${${$exonRngByGeneHsh{$geneID}}{$exonID}}{"start"} = $featureStart;
	my %refJunctStrByIntronIDHsh = %{$_[18]}; #---$junctStrByIntronIDHsh{$intronID} = $junctStr;
	my %refIntronRngByGeneHsh = %{$_[19]}; #---${${$intronRngByGeneHsh{$geneID}}{$intronID}}{"start"} = $intronStart;
	my %jClusterSSOverlapHsh = %{$_[20]}; #---${$jClusterSSOverlapHsh{$refCluster}}{$qryCluster}++ if ($refCluster ne $qryCluster);
	my %selectedIDCodonTableHsh = %{$_[21]};
	my %selectedIDStartCodonHsh = %{$_[22]};
	my %fastaHsh = %{$_[23]};
	my %nameByGeneHsh = %{$_[24]};
	my %refCtgryByGeneHsh = %{$_[25]};
	my %countFturCovInfoHsh = %{$_[26]};
	my %NGSJunctInfoHsh = %{$_[27]};
	my %SSOvrlpNGSJRefJHsh = %{$_[28]};
	my %refJunctInfoHsh = %{$_[29]};
	
	printProgressScale("Reconstructing the alternative splicing isforms based on JClusters", 50);
	
	#---go through all reference transcripts
	my $procGeneNum = my $altSplcGeneNum = my $totalIsfmNum = my $isfmCompleteORFNum = 0;
	
	open (ALTSPLCLOG, ">$outDir/altSpliceInfo/alt.splicing.allRefTrnscpt.log.txt");
	print ALTSPLCLOG join "", ((join "\t", ("refTrnscptID", "ctgry", "alt.altNum", "altFullORF", "exonSkipExactType", "skippedJunctNum", "skippeRefJNum", "altJCluster", "altSpliceTypeHsh{exonSkip}", "altSpliceTypeHsh{intronCreate}", "altSpliceTypeHsh{others}", "skippedJCluster", "altJClusterReadNum", "altJNGSConfirmedJReadNumRatio", "covPerNtRefIsfm", "senseSplicingEfficiency")), "\n");
	open (REFGFF, ">$outDir/GFFAndBED/reference.only.gff");
	open (BOTHGFF, ">$outDir/GFFAndBED/both.reference.alternative.isoform.gff");
	open (ALTGFF, ">$outDir/GFFAndBED/alternative.isoform.only.gff");
	open (ALTJUNCLOG, ">$outDir/altSpliceInfo/alternative.junction.log.txt");
	print ALTJUNCLOG "refTrnscptID.alt"."altNum"."\t"."allJClusterStr"."\n";
	print REFGFF "##gff-version\t3\n";
	print BOTHGFF "##gff-version\t3\n";
	print ALTGFF "##gff-version\t3\n";
	
	my $totalGeneNum = keys %refExonRngByGeneHsh;
	
	my (%altTrnscptInfoHsh, %altJunctStrInfoHsh, %altSplicingGeneCount, %exonNumGeneCountHsh, %exonSkipOrintronCreateGeneCountHsh, %NGSConfirmedRefJunctInfoHsh);
	
	open (INTRONCREATECOVPERNT, ">$outDir/correlation/intronCreate.gene.covPerNT.log.txt");
	open (NOTINTRONCREATECOVPERNT, ">$outDir/correlation/nonintronCreate.gene.covPerNT.log.txt");

	#---go through each cntg
	foreach my $cntg (sort {$a cmp $b} keys %refRngXSHsh) {
		
		my $cntgSeq = $fastaHsh{$cntg};
		
		#---go through each gene on cntg
		foreach my $refTrnscptID (sort {$a cmp $b} keys %{$refRngXSHsh{$cntg}}) {
			$procGeneNum++;

			my $exonLength = 0;
			foreach my $exonID (keys %{$refExonRngByGeneHsh{$refTrnscptID}}) {
				$exonLength += ${${$refExonRngByGeneHsh{$refTrnscptID}}{$exonID}}{"end"} - ${${$refExonRngByGeneHsh{$refTrnscptID}}{$exonID}}{"start"};
			}
			
			my $ctgry = $refCtgryByGeneHsh{$refTrnscptID};
			updateProgressBar($refTrnscptID, $procGeneNum, $totalGeneNum, 50, 10);
			my $strnd = $refStrndHsh{$refTrnscptID};
			my $covPerNt = 0;
			
			if (exists $countFturCovInfoHsh{$refTrnscptID}) {

				if ($strndSpecificPileup eq "yes") {
					if ($strnd eq "+") {
						$covPerNt = sprintf "%.06f", ${$countFturCovInfoHsh{$refTrnscptID}}{"plusCov"}/${$countFturCovInfoHsh{$refTrnscptID}}{"length"};
					} else {
						$covPerNt = sprintf "%.06f", ${$countFturCovInfoHsh{$refTrnscptID}}{"minusCov"}/${$countFturCovInfoHsh{$refTrnscptID}}{"length"};
					}
				} else {
					$covPerNt = sprintf "%.06f", (${$countFturCovInfoHsh{$refTrnscptID}}{"plusCov"} + ${$countFturCovInfoHsh{$refTrnscptID}}{"minusCov"})/${$countFturCovInfoHsh{$refTrnscptID}}{"length"};
				}
			}
			
			#---get all the refJunctStr of the transcript
			my @refJunctStrAry;
			my %refJInNGSJClusterHsh;
			my @tmpNGSConfirmedRefJunctReadNumAry; #---to contain all supported read nums for all ref junct

			#---NGSJCluster will be classified into one of the three categories base on their overlapping with ref
			%{$refJInNGSJClusterHsh{"ovrlpCmplt"}} = (); #--- JCluster contains a NGSJunct that is completely overlapping with the refJunct;
			%{$refJInNGSJClusterHsh{"ovrlpInferior"}} = (); #--- JCluster contains a NGSJunct that is overlapping with the refJunct but not a super junction;
			%{$refJInNGSJClusterHsh{"alternative"}} = (); #--- Any JCluster that is doesnt belong to any of the above two, which will be considered as alternative junction cluster

			my $intronCreate = "n";

			#---check if refJunct is confirmed by NGS junct
			#---if there's refJunctStr on refTrnscptID
			if (exists $refIntronRngByGeneHsh{$refTrnscptID}) {
				foreach my $refJunctStr (keys %{$refIntronRngByGeneHsh{$refTrnscptID}}) {
					push @refJunctStrAry, $refJunctStr;

					${$refJunctInfoHsh{$refJunctStr}}{"confirmedNGSJReadNum"} = "none";
					${$refJunctInfoHsh{$refJunctStr}}{"ovlpNGSJInferior"} = "none";
					${$refJunctInfoHsh{$refJunctStr}}{"refGeneCovPerNt"} = $covPerNt;
					${$refJunctInfoHsh{$refJunctStr}}{"readNumCovPerNtRatio"} = "none";

					#---if there is NGSJunctStr overlapping with refJunctStr
					if (exists $SSOvrlpRefJNGSJHsh{$refJunctStr}) {
						${$refJunctInfoHsh{$refJunctStr}}{"ovlpNGSJInferior"} = "no";

						foreach my $NGSJunctStr (keys %{$SSOvrlpRefJNGSJHsh{$refJunctStr}}) {

							#---get all info of the NGSJunct
							my @NGSJunctStrSplt = split /:/, $NGSJunctStr;
							my $intronSize = $NGSJunctStrSplt[2] - $NGSJunctStrSplt[1] + 1;
							my $jCluster = $jClusterNameByJunctHsh{$NGSJunctStr};
							my $readNum = ${$NGSJunctInfoHsh{$NGSJunctStr}}{"readNum"};
							my $superiority = "inferior";
							$superiority = "super" if (exists $superJunctOvrlpClusterHsh{$NGSJunctStr});

							#----any overlapping junction that is inferior
							${$refJunctInfoHsh{$refJunctStr}}{"ovlpNGSJInferior"} = "yes" if ($superiority eq "inferior");

							#---if NGSJunctStr and refJunctStr overlaps completely, i.e. scenerio zero
							if (${$SSOvrlpRefJNGSJHsh{$refJunctStr}}{$NGSJunctStr} == 0) {
								push @tmpNGSConfirmedRefJunctReadNumAry, $readNum;
								${$NGSConfirmedRefJunctInfoHsh{$NGSJunctStr}}{"readNum"} = $readNum; 
								${$NGSConfirmedRefJunctInfoHsh{$NGSJunctStr}}{"superiority"} = $superiority;
								${$NGSConfirmedRefJunctInfoHsh{$NGSJunctStr}}{"intronSize"} = $intronSize;
								${$refJInNGSJClusterHsh{"ovrlpCmplt"}}{$jCluster} = $superiority;
								${$refJunctInfoHsh{$refJunctStr}}{"confirmedNGSJReadNum"} = $readNum;
								${$refJunctInfoHsh{$refJunctStr}}{"readNumCovPerNtRatio"} = sprintf "%.06f", $readNum/$covPerNt;

							#---if NGSJunctStr and refJunctStr not overlaps completely but is inferior
							} elsif ($superiority eq "inferior") {
								${$refJInNGSJClusterHsh{"ovrlpInferior"}}{$jCluster} = $superiority;
							}
						}
					}
				}
			}

			#---calculated the mean and SD of read num of the NGS confirmed refJunctions
			my $NGSConfirmedRefJReadNumMean = my $NGSConfirmedRefJReadNumSD = "null";
			my $NGSConfirmedRefJNum = @tmpNGSConfirmedRefJunctReadNumAry;
			${$exonNumGeneCountHsh{$NGSConfirmedRefJNum}}{$refTrnscptID}++;
			if (@tmpNGSConfirmedRefJunctReadNumAry > 0) {
				($NGSConfirmedRefJReadNumMean, $NGSConfirmedRefJReadNumSD) = calculateStandardDeviationAndMean(\@tmpNGSConfirmedRefJunctReadNumAry);
			}

			my $refIntronNum = @refJunctStrAry;
			my $NGSConfirmedJunctNum = keys %{$refJInNGSJClusterHsh{"ovrlpCmplt"}};

			#---get all the junction clusters on the reference transcripts
			my %jClusterOnRefHsh;#---define empty hashes for the convenience of counting later on
			foreach my $senseness (("sense", "antisense")) {
				foreach my $superiority (("super", "inferior")) {
					%{${${$jClusterOnRefHsh{$refTrnscptID}}{$senseness}}{$superiority}} = (); #---define the hash
				}
			}

			my %altJClusterHsh;
			
			##############cluster based reconstruction############################
			#---If the refTranscript is overlapping with NGSJunct
			if (exists $XSOvrlpRefTrnscptNGSJHsh{$refTrnscptID}) {

				#---goes through all NGSJunct
				foreach my $NGSJunctStr (keys %{$XSOvrlpRefTrnscptNGSJHsh{$refTrnscptID}}) {

					#-----record the junction only if the whole junct range is within the refTrnscpt rng
					if (${$XSOvrlpRefTrnscptNGSJHsh{$refTrnscptID}}{$NGSJunctStr} == 3) {

						my $jCluster = $jClusterNameByJunctHsh{$NGSJunctStr};
						my $jClusterProminentJStr = ${$jClusterInfoHsh{$jCluster}}{"prominentJunct"};
						
						#----make sure the prominent junct is inside the gene model
						if (${$XSOvrlpRefTrnscptNGSJHsh{$refTrnscptID}}{$jClusterProminentJStr} == 3) {
							my $superiority = "inferior";
							$superiority = "super" if (exists $superJunctOvrlpClusterHsh{$NGSJunctStr});

							#---hit in the SS hash, sense hit
							if (exists ${$SSOvrlpRefTrnscptNGSJHsh{$refTrnscptID}}{$NGSJunctStr}) {
								push @{${${${$jClusterOnRefHsh{$refTrnscptID}}{"sense"}}{$superiority}}{$jCluster}}, $NGSJunctStr;

								#---the jCluster is no complete overlap (ovrlpCmplt) or partial inferior overlap (ovrlpInferior);
								if ((not exists ${$refJInNGSJClusterHsh{"ovrlpCmplt"}}{$jCluster}) and (not exists ${$refJInNGSJClusterHsh{"ovrlpInferior"}}{$jCluster})) {
									${$refJInNGSJClusterHsh{"alternative"}}{$jCluster} = $superiority;
								}

							#---antisense hit
							} else {
								push @{${${${$jClusterOnRefHsh{$refTrnscptID}}{"antisense"}}{$superiority}}{$jCluster}}, $NGSJunctStr;
							}
						}
					}
				}
			}

			my $altJCluterNum = keys %{$refJInNGSJClusterHsh{"alternative"}};

			#print "$refTrnscptID, intron:$refIntronNum/$NGSConfirmedJunctNum   altJCluster:$altJCluterNum         \n";
			
			#---reconstruct the alternative splicing isoforms
			#
			# NOTE:!!!!!! RATIONALE!!!!!
			# Each alternative isoform is a variation of the reference isoform, with ONE single additional alternative JCluster, any prominent JCluster
			# overlapping with the alternative JCluster will be removed.
			
			#---store all exon bounds of the reference isoform
			my %allIsfmExonBoundHsh;
			my @refIsfmExonBoundAry;
			
			foreach my $exonID (keys %{$refExonRngByGeneHsh{$refTrnscptID}}) {
				push @refIsfmExonBoundAry, ${${$refExonRngByGeneHsh{$refTrnscptID}}{$exonID}}{"start"};
				push @refIsfmExonBoundAry, ${${$refExonRngByGeneHsh{$refTrnscptID}}{$exonID}}{"end"};
			}
			
			@refIsfmExonBoundAry = sort {$a<=>$b} @refIsfmExonBoundAry;
			my $geneStart = $refIsfmExonBoundAry[0];
			my $geneEnd = $refIsfmExonBoundAry[-1];
			@{${$allIsfmExonBoundHsh{$refTrnscptID}}{"ref"}} = @refIsfmExonBoundAry;
			
			my $totalAltIsfmNum = keys %{$refJInNGSJClusterHsh{"alternative"}};

			#---get the seq of the ref
			my ($refDNASeq, $refFirst3Bases, $refLast3Bases, $refPepSeq, $refFullORF) = extractExonBoundSeqAndTranslate(\@refIsfmExonBoundAry, $cntgSeq, $strnd, \%selectedIDCodonTableHsh);
			
			#---generate GFF files
			my $refGenemRNALineStr = generatemRNAExonGFFLine($refTrnscptID, $refTrnscptID.".ref", \@refIsfmExonBoundAry, $cntg, $strnd, "y", $nameByGeneHsh{$refTrnscptID}, "completeORF=".$refFullORF, "y");
			my $refGeneOnlyLineStr = generatemRNAExonGFFLine($refTrnscptID, $refTrnscptID.".ref", \@refIsfmExonBoundAry, $cntg, $strnd, "y", $nameByGeneHsh{$refTrnscptID}, "completeORF=".$refFullORF, "n");

			print REFGFF $refGenemRNALineStr."\n";
			print BOTHGFF $refGenemRNALineStr."\n";
			print ALTGFF $refGeneOnlyLineStr."\n" if ($totalAltIsfmNum > 0);
			
			my $altNum = 0;
			
			my %altSpliceTypeHsh;
			
			#---goes through all alternative jClusters
			foreach my $altJCluster (keys %{$refJInNGSJClusterHsh{"alternative"}}) {
				
				$altSpliceTypeHsh{"exonSkip"} = "n"; #----the alt junct is a super junct but it's ovrlapping junction is not a prominent junction
				$altSpliceTypeHsh{"intronCreate"} = "n";	#----intron creation
				$altSpliceTypeHsh{"others"} = "n";	#----other scinerios
			
				$altNum++;
				$totalIsfmNum++;
				my $altJClusterReadNum = ${$jClusterInfoHsh{$altJCluster}}{"mostSupportReadNum"};
				my $altJNGSConfirmedJReadNumRatio = "null";
				$altJNGSConfirmedJReadNumRatio = sprintf "%.06f", $altJClusterReadNum/$NGSConfirmedRefJReadNumMean if ($NGSConfirmedRefJReadNumMean ne "null");
				
				my $altJunctStr = ${$jClusterInfoHsh{$altJCluster}}{"prominentJunct"};
				my $exonSkipExactType = "null";
				$exonSkipExactType = $totalExonSkippingTypeHsh{$altJunctStr} if (exists $totalExonSkippingTypeHsh{$altJunctStr});
				my %skippeRefJunctHsh;
				my @skippeRefJunctAry;
				my $skippeRefJNum = 0;
				
				#---store refJunctStr and the altJunctStr in a hash
				my %allJunctStrOnIsoformHsh;
				$allJunctStrOnIsoformHsh{$altJunctStr}++;
				
				foreach my $refJunctStr (@refJunctStrAry) {
					$allJunctStrOnIsoformHsh{$refJunctStr}++;
				}
				
				#---check for overlapping with the prominent JClusters, remove pJCluster if it overlaps with altJCluster
				#---SKIPPING THE JUNCTIONS HERE!

				#---if the alt junct is overlapping with some junct
				if (exists $SSOvrlpNGSJRefJHsh{$altJunctStr}) {
					
					#---go through all overlapping jCluster
					foreach my $ovrlapRefJunct (keys %{$SSOvrlpNGSJRefJHsh{$altJunctStr}}) {
						
						#---if the alt junct is a super junct
						if (exists $superJunctOvrlpClusterHsh{$altJunctStr}) {
							$skippeRefJNum++;
							$altSpliceTypeHsh{"exonSkip"} = "y";
							$skippeRefJunctHsh{$ovrlapRefJunct}++;
							push @skippeRefJunctAry, $ovrlapRefJunct;
							delete $allJunctStrOnIsoformHsh{$ovrlapRefJunct}; #---remove of the overlapped ref junction
						} else {
							$altSpliceTypeHsh{"others"} = "y";
						}
					}

				} else {#--not overlap with refJunct
					$altSpliceTypeHsh{"intronCreate"} = "y";
				}
				
				my $skippedJunctNum = keys %skippeRefJunctHsh;
				push @skippeRefJunctAry, "na" if (@skippeRefJunctAry == 0);
				my $skippedJCluster = join ";", @skippeRefJunctAry;
				
				#---get all the intron bounds and store them
				my @altIsfmExonBoundAry;
				push @altIsfmExonBoundAry, ($geneStart, $geneEnd);
				my $allJunctStr = "";
				foreach my $junctStr (keys %allJunctStrOnIsoformHsh) {
					my @junctStrSplt = split /:/, $junctStr;
					my $junctBoundStart = $junctStrSplt[1]-1;
					my $junctBoundEnd = $junctStrSplt[2]+1;
					push @altIsfmExonBoundAry, ($junctBoundStart, $junctBoundEnd);
					$allJunctStr .= $junctStr."\t";
				}
				
				@altIsfmExonBoundAry = sort {$a<=>$b} @altIsfmExonBoundAry;
				@{${$allIsfmExonBoundHsh{$refTrnscptID}}{"alt".$altNum}} = @altIsfmExonBoundAry;
				
				print ALTJUNCLOG $refTrnscptID.".alt".$altNum."\t".$allJunctStr."\n";
				my ($altDNASeq, $altFirst3Bases, $altLast3Bases, $altPepSeq, $altFullORF) = extractExonBoundSeqAndTranslate(\@altIsfmExonBoundAry, $cntgSeq, $strnd, \%selectedIDCodonTableHsh);
				my $altmRNAOnlyLineStr = generatemRNAExonGFFLine($refTrnscptID, $refTrnscptID.".alt".$altNum, \@altIsfmExonBoundAry, $cntg, $strnd, "n", $nameByGeneHsh{$refTrnscptID}, "completeORF=".$altFullORF, "y");
				print BOTHGFF $altmRNAOnlyLineStr."\n";
				print ALTGFF $altmRNAOnlyLineStr."\n";
				$isfmCompleteORFNum++ if ($altFullORF eq "y");
				
				${${$altTrnscptInfoHsh{$refTrnscptID}}{$altNum}}{"altFullORF"} = $altFullORF;
				${${$altTrnscptInfoHsh{$refTrnscptID}}{$altNum}}{"skippedJunctNum"} = $skippedJunctNum;
				${${$altTrnscptInfoHsh{$refTrnscptID}}{$altNum}}{"skippeRefJNum"} = $skippeRefJNum;
				${${$altTrnscptInfoHsh{$refTrnscptID}}{$altNum}}{"altJCluster"} = $altJCluster;
				${${$altTrnscptInfoHsh{$refTrnscptID}}{$altNum}}{"exonSkip"} = $altSpliceTypeHsh{"exonSkip"};
				${${$altTrnscptInfoHsh{$refTrnscptID}}{$altNum}}{"intronCreate"} = $altSpliceTypeHsh{"intronCreate"};
				${${$altTrnscptInfoHsh{$refTrnscptID}}{$altNum}}{"otherAltType"} = $altSpliceTypeHsh{"others"};
				${${$altTrnscptInfoHsh{$refTrnscptID}}{$altNum}}{"skippedJCluster"} = $skippedJCluster;
				${${$altTrnscptInfoHsh{$refTrnscptID}}{$altNum}}{"altJClusterReadNum"} = $altJClusterReadNum;
				${${$altTrnscptInfoHsh{$refTrnscptID}}{$altNum}}{"altJNGSConfirmedJReadNumRatio"} = $altJNGSConfirmedJReadNumRatio;
				${${$altTrnscptInfoHsh{$refTrnscptID}}{$altNum}}{"ctgry"} = $ctgry;
				${${$altTrnscptInfoHsh{$refTrnscptID}}{$altNum}}{"exonSkipExactType"} = $exonSkipExactType;
				${${$altTrnscptInfoHsh{$refTrnscptID}}{$altNum}}{"altJunctStr"} = $altJunctStr;
				my $senseSplicingEfficiency = ${$NGSJunctInfoHsh{$altJunctStr}}{"senseSplicingEfficiency"};
				${${$altTrnscptInfoHsh{$refTrnscptID}}{$altNum}}{"senseSplicingEfficiency"} = $senseSplicingEfficiency;
				
				${$altJunctStrInfoHsh{$altJunctStr}}{"exonSkipExactType"} = $exonSkipExactType;
				${$altJunctStrInfoHsh{$altJunctStr}}{"otherAltType"} = $altSpliceTypeHsh{"others"};
				${$altJunctStrInfoHsh{$altJunctStr}}{"intronCreate"} = $altSpliceTypeHsh{"intronCreate"};
				${$altJunctStrInfoHsh{$altJunctStr}}{"exonSkip"} = $altSpliceTypeHsh{"exonSkip"};

				if ($skippeRefJNum >= 2) {

					${$exonSkipOrintronCreateGeneCountHsh{"exonSkip"}}{$refTrnscptID}++;

					if ($exonSkipExactType eq "exact") {
						${$altSplicingGeneCount{"exactSkipTwoPJMore"}}{$refTrnscptID}++;
					} else {
						${$altSplicingGeneCount{"inexactSkipTwoPJMore"}}{$refTrnscptID}++;
					}
				} 
				
				if (($skippedJunctNum >= 1) and ($skippeRefJNum < 2)) {	
					${$exonSkipOrintronCreateGeneCountHsh{"intronCreate"}}{$refTrnscptID}++;

					if ($exonSkipExactType eq "exact") {
						${$altSplicingGeneCount{"exactSkipTwoPJLess"}}{$refTrnscptID}++;
					} else {
						${$altSplicingGeneCount{"inexactSkipTwoPJLess"}}{$refTrnscptID}++;
					}
				}
				
				if ($altSpliceTypeHsh{"others"} eq "y") {
					${$altSplicingGeneCount{"otherAltType"}}{$refTrnscptID}++;
					${$exonSkipOrintronCreateGeneCountHsh{"intronCreate"}}{$refTrnscptID}++;
				}
			
				$intronCreate = "y" if (($altSpliceTypeHsh{"intronCreate"} eq "y") and ($altJClusterReadNum >= 1));#----adhoc cuttoff
				
				print ALTSPLCLOG join "", ((join "\t", ($refTrnscptID, $ctgry, "alt".$altNum, $altFullORF, $exonSkipExactType, $skippedJunctNum, $skippeRefJNum, $altJCluster, $altSpliceTypeHsh{"exonSkip"}, $altSpliceTypeHsh{"intronCreate"}, $altSpliceTypeHsh{"others"}, $skippedJCluster, $altJClusterReadNum, $altJNGSConfirmedJReadNumRatio, $covPerNt, $senseSplicingEfficiency)), "\n");
			}
			
			if ($intronCreate eq "y") {
				print INTRONCREATECOVPERNT $covPerNt."\n";
			} else {
				print NOTINTRONCREATECOVPERNT $covPerNt."\n";
			}
			
			$altSplcGeneNum++ if ($altNum > 0);
		}
	}
	
	print "\n\n";
	close ALTSPLCLOG;
	close BOTHGFF;
	close REFGFF;
	close ALTGFF;
	
	$totalIsfmNum += $altSplcGeneNum;

	#---calculate constitutive junction size and readnnm stats
	my (@tmpCnstvJunctReadNum, @tmpCnstvJunctSize);
	foreach my $NGSJunctStr (keys %NGSConfirmedRefJunctInfoHsh) {
		push @tmpCnstvJunctReadNum, ${$NGSConfirmedRefJunctInfoHsh{$NGSJunctStr}}{"readNum"}; 
		push @tmpCnstvJunctSize, ${$NGSConfirmedRefJunctInfoHsh{$NGSJunctStr}}{"intronSize"}; 
	}
	my $cnstivJunctNum = keys %NGSConfirmedRefJunctInfoHsh;
	my ($allCnstivJReadNumMean, $allCnstivJReadNumSD) = calculateStandardDeviationAndMean(\@tmpCnstvJunctReadNum);
	my ($allCnstivJIntronSizeMean, $allCnstivJIntronSizeSD) = calculateStandardDeviationAndMean(\@tmpCnstvJunctSize);
	
	my $isfmPerAltSplcGene = sprintf "%.06f", $totalIsfmNum/$altSplcGeneNum;
	my $exactSkipTwoPJMore = keys %{$altSplicingGeneCount{"exactSkipTwoPJMore"}};
	my $inexactSkipTwoPJMore = keys %{$altSplicingGeneCount{"inexactSkipTwoPJMore"}};
	my $exactSkipTwoPJLess = keys %{$altSplicingGeneCount{"exactSkipTwoPJLess"}};
	my $inexactSkipTwoPJLess = keys %{$altSplicingGeneCount{"inexactSkipTwoPJLess"}};
	my $otherAltType = keys %{$altSplicingGeneCount{"otherAltType"}};
	my $exonSkipNum = keys %{$exonSkipOrintronCreateGeneCountHsh{"exonSkip"}};
	my $intronCreateNum = keys %{$exonSkipOrintronCreateGeneCountHsh{"intronCreate"}};

	print "Total number of genes = $procGeneNum\n";
	print "Number of consitutive junctions = $cnstivJunctNum\n";
	print "Mean and SD consitutive junctions Read Num = $allCnstivJReadNumMean, $allCnstivJReadNumSD\n";
	print "Mean and SD consitutive junctions Intron Size = $allCnstivJIntronSizeMean, $allCnstivJIntronSizeSD\n";
	print "Number of alt spliced genes = $altSplcGeneNum\n";
	foreach my $NGSConfirmedRefJNum (sort {$a <=> $b} keys %exonNumGeneCountHsh) {
		my $tmpGeneNum = keys %{$exonNumGeneCountHsh{$NGSConfirmedRefJNum}};
		print "Number of genes with $NGSConfirmedRefJNum constitutive junctions = $tmpGeneNum\n";
	}
	
	my %exonSkipOrintronCreateGeneOverlapHsh;
	open (CONSTITIVEXONSKIP, ">$outDir/altSpliceInfo/exonSkip.genelist.txt");
	foreach my $geneID (keys %{$exonSkipOrintronCreateGeneCountHsh{"exonSkip"}}) {
		my $tmpOvrlp;
		if (exists ${$exonSkipOrintronCreateGeneCountHsh{"intronCreate"}}{$geneID}) {
			$tmpOvrlp = "both";
		} else {
			$tmpOvrlp = "exonSkipOnly";
		}
		${$exonSkipOrintronCreateGeneOverlapHsh{$tmpOvrlp}}{$geneID}++;
		print CONSTITIVEXONSKIP $geneID."\t".$tmpOvrlp."\n";
	}
	
	open (CONSTITIVINTRONRETN, ">$outDir/altSpliceInfo/intronCreate.genelist.txt");
	foreach my $geneID (keys %{$exonSkipOrintronCreateGeneCountHsh{"intronCreate"}}) {
		my $tmpOvrlp;
		if (exists ${$exonSkipOrintronCreateGeneCountHsh{"exonSkip"}}{$geneID}) {
			$tmpOvrlp = "both";
		} else {
			$tmpOvrlp = "intronCreateOnly";
		}
		${$exonSkipOrintronCreateGeneOverlapHsh{$tmpOvrlp}}{$geneID}++;
		print CONSTITIVINTRONRETN $geneID."\t".$tmpOvrlp."\n";
	}

	my $intronCreateOnly = keys %{$exonSkipOrintronCreateGeneOverlapHsh{"intronCreateOnly"}};
	my $bothintronCreateExonSkip = keys %{$exonSkipOrintronCreateGeneOverlapHsh{"both"}};
	my $exonSkipOnly = keys %{$exonSkipOrintronCreateGeneOverlapHsh{"exonSkipOnly"}};

	print "Number of isoform per alt spliced gene = $isfmPerAltSplcGene\n";
	print "Total number of alternative isoforms = $totalIsfmNum\n";
	print "Number of alternative isoforms with complete ORF= $isfmCompleteORFNum\n";
	print "Number of genes with exact skipping on > 2 reference junctions= $exactSkipTwoPJMore\n";
	print "Number of genes with inexact skipping on > 2 reference junctions= $inexactSkipTwoPJMore\n";
	print "Number of genes with exact skipping on < 2 reference junctions= $exactSkipTwoPJLess\n";
	print "Number of genes with inexact skipping on < 2 reference junctions= $inexactSkipTwoPJLess\n";
	print "Number of genes with non-exon-skipping type of alt splicing= $otherAltType\n";
	print "Number of genes with exonSkip = $exonSkipNum\n";
	print "Number of genes with intronCreate= $intronCreateNum\n";
	print "Number of genes with both exonSkip and intronCreate = $bothintronCreateExonSkip\n";
	print "Number of genes with only exonSkip = $exonSkipOnly\n";
	print "Number of genes with only intronCreate = $intronCreateOnly\n";
	
	printProgressScale("Reconstructing the alternative splicing isforms based on individual junctions", 50);
	
	#---go through all reference transcripts
	my $indivJProcGeneNum = my $indivJAltSplcGeneNum = my $indivJTotalIsfmNum = my $indivJIsfmCompleteORFNum = 0;

	#######################Junction based reconstruction, ignore the clusters, each junction can be an alternative junction

	open (INDIVJREFGFF, ">$outDir/GFFAndBED/indivJ.reference.only.gff");
	open (INDIVJBOTHGFF, ">$outDir/GFFAndBED/indivJ.both.reference.alternative.isoform.gff");
	open (INDIVJALTGFF, ">$outDir/GFFAndBED/indivJ.alternative.isoform.only.gff");
	print INDIVJREFGFF "##gff-version\t3\n";
	print INDIVJBOTHGFF "##gff-version\t3\n";
	print INDIVJALTGFF "##gff-version\t3\n";

	my %junctBasedAltIsofmInfoHsh;

	#---set the default values of the alternative splicing info
	foreach my $NGSJunctStr (keys %NGSJunctInfoHsh) {
		${$NGSJunctInfoHsh{$NGSJunctStr}}{"refOrAlt"} = "none"; #----none, alt or ref;
		${$NGSJunctInfoHsh{$NGSJunctStr}}{"isofmID"} = "none"; #----none, isofmID or ref;
		${$NGSJunctInfoHsh{$NGSJunctStr}}{"altIsofmFullORF"} = "none"; #----none, ref, yes or no to preserved the original ORF of the refIsofm;
		${$NGSJunctInfoHsh{$NGSJunctStr}}{"altSplcType"} = "none"; #----none, ref,intronCreate, altSite or exonSkip
		${$NGSJunctInfoHsh{$NGSJunctStr}}{"onRefGene"} = "none"; #----none, or gene name;
		${$NGSJunctInfoHsh{$NGSJunctStr}}{"refJunctAvgReadNum"} = "none"; #----none, mean number of reads of the reference Junctions
		${$NGSJunctInfoHsh{$NGSJunctStr}}{"NGSRefJunctReadNumRatio"} = "none"; #----none, ratio of NGS junct readNum to the meanRefJunctReadNum
		${$NGSJunctInfoHsh{$NGSJunctStr}}{"readNumCovPerNtRatio"} = "none"; #----none, or readNum to CovPerNt of the refGene Ratio
		${$NGSJunctInfoHsh{$NGSJunctStr}}{"altSiteType"} = "none"; #----none, or ref, altSite5 or altSite3 or altSiteBoth
		${$NGSJunctInfoHsh{$NGSJunctStr}}{"exonSkipType"} = "none"; #----none, or ref, incomplete or complete or partial
		${$NGSJunctInfoHsh{$NGSJunctStr}}{"refShift5"} = "none"; #----none, or ref, alternative site/exon Skip shift from ref J at 5 end
		${$NGSJunctInfoHsh{$NGSJunctStr}}{"refShift3"} = "none"; #----none, or ref, alternative site/exon Skip shift from ref J at 3 end
		${$NGSJunctInfoHsh{$NGSJunctStr}}{"skippedExonNum"} = "none"; #----none, or ref, number of skipped exons
		${$NGSJunctInfoHsh{$NGSJunctStr}}{"ovrlpRefJunctRdRatio"} = "none"; #----none, ref, or ratio of the superJunction rd num versus the mean readNum of skipped the inferior Junctions
		${$NGSJunctInfoHsh{$NGSJunctStr}}{"avgOvrlpRefJReadNum"} = "none"; #----none, ref, or average readNum of the overlapping reference junction
	}

	#---go through each cntg
	open (FASTASEQ, ">$outDir/indivJunctAltIsfmSeq.fasta");
	foreach my $cntg (sort {$a cmp $b} keys %refRngXSHsh) {
		
		my $cntgSeq = $fastaHsh{$cntg};
		
		#---go through each gene on cntg
		foreach my $refTrnscptID (sort {$a cmp $b} keys %{$refRngXSHsh{$cntg}}) {
			$indivJProcGeneNum++;

			my $ctgry = $refCtgryByGeneHsh{$refTrnscptID};
			updateProgressBar($refTrnscptID, $indivJProcGeneNum, $totalGeneNum, 50, 10);
			my $strnd = $refStrndHsh{$refTrnscptID};


			my $exonLength = 0;
			foreach my $exonID (keys %{$refExonRngByGeneHsh{$refTrnscptID}}) {
				$exonLength += ${${$refExonRngByGeneHsh{$refTrnscptID}}{$exonID}}{"end"} - ${${$refExonRngByGeneHsh{$refTrnscptID}}{$exonID}}{"start"};
			}

			my $covPerNt = 0;
			
			if (exists $countFturCovInfoHsh{$refTrnscptID}) {

				if ($strndSpecificPileup eq "yes") {
					if ($strnd eq "+") {
						$covPerNt = sprintf "%.06f", ${$countFturCovInfoHsh{$refTrnscptID}}{"plusCov"}/${$countFturCovInfoHsh{$refTrnscptID}}{"length"};
					} else {
						$covPerNt = sprintf "%.06f", ${$countFturCovInfoHsh{$refTrnscptID}}{"minusCov"}/${$countFturCovInfoHsh{$refTrnscptID}}{"length"};
					}
				} else {
					$covPerNt = sprintf "%.06f", (${$countFturCovInfoHsh{$refTrnscptID}}{"plusCov"} + ${$countFturCovInfoHsh{$refTrnscptID}}{"minusCov"})/${$countFturCovInfoHsh{$refTrnscptID}}{"length"};
				}
			}

			
			my %altIsfmJunctStrHsh;
			my $isofmNum = 0;
			
			#---get the NGS confirmed RefJunct info
			#my %NGSConirmedRefJInfoHsh;
			my $confirmedRefJNum = 0;
			my $confirmedRefJReadNum = 0;
			my $refJunctNum = keys %{$refIntronRngByGeneHsh{$refTrnscptID}};
			foreach my $refJunctStr (keys %{$refIntronRngByGeneHsh{$refTrnscptID}}) {
				#---refJ Confirmed by NGSJ
				if (exists $NGSJunctInfoHsh{$refJunctStr}) {
					$confirmedRefJNum++;
					my $readNum = ${$NGSJunctInfoHsh{$refJunctStr}}{"readNum"};
					$confirmedRefJReadNum += $readNum;
				}
			}
			
			my $refJunctAvgReadNum = "none";
			if ($refJunctNum > 0) {
				$refJunctAvgReadNum = 0;
				if ($confirmedRefJReadNum > 0) {
					$refJunctAvgReadNum = sprintf "%.06f", $confirmedRefJReadNum/$confirmedRefJNum;
				}
			}
			
			#---If the refTranscript is overlapping with NGSJunct
			if (exists $XSOvrlpRefTrnscptNGSJHsh{$refTrnscptID}) {

				#---goes through all NGSJunct, each NGSJunct is an alternative isoform
				foreach my $NGSJunctStr (keys %{$XSOvrlpRefTrnscptNGSJHsh{$refTrnscptID}}) {
					
					#-----record the junction only if the whole junct range is within the refTrnscpt rng
					if (${$XSOvrlpRefTrnscptNGSJHsh{$refTrnscptID}}{$NGSJunctStr} == 3) {
						
						#---hit in the SS hash, sense hit
						if (exists ${$SSOvrlpRefTrnscptNGSJHsh{$refTrnscptID}}{$NGSJunctStr}) {

							my @NGSJunctStrSplt = split /:/, $NGSJunctStr;
							my $NGSJunctBoundStart = $NGSJunctStrSplt[1]-1;
							my $NGSJunctBoundEnd = $NGSJunctStrSplt[2]+1;

							my $readNum = ${$NGSJunctInfoHsh{$NGSJunctStr}}{"readNum"};
							my $NGSRefJunctReadNumRatio = "none";
							if ($confirmedRefJReadNum > 0) {
								$NGSRefJunctReadNumRatio = sprintf "%.06f", $readNum/$refJunctAvgReadNum;
							}
							my $readNumCovPerNtRatio =  sprintf "%.06f", $readNum/$covPerNt;

							${$NGSJunctInfoHsh{$NGSJunctStr}}{"onRefGene"} = $refTrnscptID; #----none, or gene name;
							${$NGSJunctInfoHsh{$NGSJunctStr}}{"refJunctAvgReadNum"} = $refJunctAvgReadNum; #----none, mean number of reads of the reference Junctions
							${$NGSJunctInfoHsh{$NGSJunctStr}}{"NGSRefJunctReadNumRatio"} = $NGSRefJunctReadNumRatio; #----none, ratio of NGS junct readNum to the refJunctAvgReadNum
							${$NGSJunctInfoHsh{$NGSJunctStr}}{"readNumCovPerNtRatio"} = $readNumCovPerNtRatio; #----none, or readNum to CovPerNt of the refGene Ratio

							#---exact hit of the refJ
							if (exists ${$refIntronRngByGeneHsh{$refTrnscptID}}{$NGSJunctStr}) {
								${$NGSJunctInfoHsh{$NGSJunctStr}}{"refOrAlt"} = "ref"; #----none, alt or ref;
								${$NGSJunctInfoHsh{$NGSJunctStr}}{"isofmID"} = "ref"; #----none, isofmID or ref;
								${$NGSJunctInfoHsh{$NGSJunctStr}}{"altIsofmFullORF"} = "ref"; #----none, ref, yes or no to preserved the original ORF of the refIsofm;
								${$NGSJunctInfoHsh{$NGSJunctStr}}{"altSplcType"} = "ref"; #----none, ref,intronCreate, altSite or exonSkip
								${$NGSJunctInfoHsh{$NGSJunctStr}}{"ovrlpRefJunctRdRatio"} = "ref"; #----none, ref, or ratio of the superJunction rd num versus the mean readNum of skipped the inferior Junctions
								${$NGSJunctInfoHsh{$NGSJunctStr}}{"avgOvrlpRefJReadNum"} = "ref"; #----none, ref, or average readNum of the overlapping reference junction
								${$NGSJunctInfoHsh{$NGSJunctStr}}{"refShift5"} = "ref"; #----none, or ref, alternative site shift at 5 end
								${$NGSJunctInfoHsh{$NGSJunctStr}}{"refShift3"} = "ref"; #----none, or ref, alternative site shift at 3 end

								${$NGSJunctInfoHsh{$NGSJunctStr}}{"altSiteType"} = "ref"; #----none, or ref, altSite5 or altSite3 or altSiteBoth

								${$NGSJunctInfoHsh{$NGSJunctStr}}{"exonSkipType"} = "ref"; #----none, or ref, incomplete or complete or partial
								${$NGSJunctInfoHsh{$NGSJunctStr}}{"skippedExonNum"} = "ref"; #----none, or ref, number of skipped exons


							#---altJ, will be on alt isfm
							} else {

								$isofmNum++;
								my $isofmID = "alt".$isofmNum;

								${$NGSJunctInfoHsh{$NGSJunctStr}}{"refOrAlt"} = "alt"; #----none, alt or ref;
								${$NGSJunctInfoHsh{$NGSJunctStr}}{"isofmID"} = $isofmID; #----none, isofmID or ref;

								${${$junctBasedAltIsofmInfoHsh{$refTrnscptID}}{$isofmID}}{"NGSJunctStr"} = $NGSJunctStr;

								#---colleact all refJunct first
								foreach my $refJunctStr (keys %{$refIntronRngByGeneHsh{$refTrnscptID}}) {
									${$altIsfmJunctStrHsh{$isofmID}}{$refJunctStr}++;
								}
								
								#---add the NGSJunct also
								${$altIsfmJunctStrHsh{$isofmID}}{$NGSJunctStr}++;
	
								#---if there is NGSJunctStr overlapping with refJunctStr, maybe exonSkip or altSplicing Site
								if (exists $SSOvrlpNGSJRefJHsh{$NGSJunctStr}) {
								
									my $ovrlapRefJunctNum = keys %{$SSOvrlpNGSJRefJHsh{$NGSJunctStr}};
	
									#---remove all overlapped refJunct from isoform
									my $allOvrlpRefJReadNum = 0;
									my $ovrlpRefJunctRdRatio = "none";
									
									my @ovrlpRefJBoundAry; 
									foreach my $refJunctStr (keys %{$SSOvrlpNGSJRefJHsh{$NGSJunctStr}}) {
										my @refJunctStrSplt = split /:/, $refJunctStr;
										my $refJunctBoundStart = $refJunctStrSplt[1]-1;
										my $refJunctBoundEnd = $refJunctStrSplt[2]+1;
										push @ovrlpRefJBoundAry, ($refJunctBoundStart, $refJunctBoundEnd);
										
										delete ${$altIsfmJunctStrHsh{$isofmID}}{$refJunctStr};
										
										if (exists $NGSJunctInfoHsh{$refJunctStr}) {
											my $tmpOvrlpRefJReadNum = ${$NGSJunctInfoHsh{$refJunctStr}}{"readNum"};
											$allOvrlpRefJReadNum += $tmpOvrlpRefJReadNum;
										}
									}

									my $avgOvrlpRefJReadNum = sprintf "%.06f", $allOvrlpRefJReadNum/$ovrlapRefJunctNum;
									$ovrlpRefJunctRdRatio = sprintf "%.06f", $readNum/$avgOvrlpRefJReadNum if ($avgOvrlpRefJReadNum > 0);

									#---get the 3 and 5 bound of the overlapping ref junction 
									@ovrlpRefJBoundAry = sort {$a <=> $b} @ovrlpRefJBoundAry;
									my $ovrlpRefJBoundStart = $ovrlpRefJBoundAry[0];
									my $ovrlpRefJBoundEnd = $ovrlpRefJBoundAry[-1];
									
									my ($refShift5, $refShift3);
									if ($strnd eq "+") {
										$refShift5 = $NGSJunctBoundStart - $ovrlpRefJBoundStart;
										$refShift3 = $NGSJunctBoundEnd - $ovrlpRefJBoundEnd;
									} elsif ($strnd eq "-") {
										$refShift3 = -1*($NGSJunctBoundStart - $ovrlpRefJBoundStart);
										$refShift5 = -1*($NGSJunctBoundEnd - $ovrlpRefJBoundEnd);
									}
									
									${$NGSJunctInfoHsh{$NGSJunctStr}}{"refShift5"} = $refShift5; #----none, or ref, alternative site shift at 5 end
									${$NGSJunctInfoHsh{$NGSJunctStr}}{"refShift3"} = $refShift3; #----none, or ref, alternative site shift at 3 end

									#---store the ovrlpRefJunctRdRatio
									${$NGSJunctInfoHsh{$NGSJunctStr}}{"ovrlpRefJunctRdRatio"} = $ovrlpRefJunctRdRatio; #----none, ref, or ratio of the superJunction rd num versus the mean readNum of skipped the inferior Junctions / alternative

									#---exonSkip
									if ($ovrlapRefJunctNum > 1) {

										${$NGSJunctInfoHsh{$NGSJunctStr}}{"altSplcType"} = "exonSkip"; #----none, ref,intronCreate, altSite or exonSkip
										${$NGSJunctInfoHsh{$NGSJunctStr}}{"skippedExonNum"} = $ovrlapRefJunctNum; #----none, or ref, number of skipped exons

										my $exonSkipType;
										if (($refShift5 == 0) and ($refShift3 == 0)) {
											$exonSkipType = "complete";
										} elsif (($refShift5 == 0) or ($refShift3 == 0)) {
											$exonSkipType = "partial";
										} else {
											$exonSkipType = "incomplete";
										}
										
										${$NGSJunctInfoHsh{$NGSJunctStr}}{"exonSkipType"} = $exonSkipType; #----none, or ref, incomplete or complete or partial

									#---altSite
									} elsif  ($ovrlapRefJunctNum == 1) {

										${$NGSJunctInfoHsh{$NGSJunctStr}}{"altSplcType"} = "altSite"; #----none, ref,intronCreate, altSite or exonSkip

										my $altSiteType;
										if (($refShift5 != 0) and ($refShift3 != 0)) {
											$altSiteType = "altSiteBoth";
										} elsif (($refShift5 != 0) and ($refShift3 == 0)) {
											$altSiteType = "altSite5";
										} elsif (($refShift5 == 0) and ($refShift3 != 0)) {
											$altSiteType = "altSite3";
										} else {
											$altSiteType = "other";
											#print "refShift5 = $refShift5 and refShift3 = $refShift3\n";
										}

										${$NGSJunctInfoHsh{$NGSJunctStr}}{"altSiteType"} = $altSiteType; #----none, or ref, altSite5 or altSite3 or altSiteBoth

									} else {
										die "debug: Impossible scenerio2\n";
									}

								#---intronCreate
								} else {
									
									${$NGSJunctInfoHsh{$NGSJunctStr}}{"altSplcType"} = "intronCreate"; #----none, ref,intronCreate, altSite or exonSkip
									
								}
							}#---end of if ${$refIntronRngByGeneHsh{$refTrnscptID}}{$NGSJunctStr} {
						}
					}
				}
			}


			#---output the reference isoform first
			#---get the seq of the ref
			my @refIsfmExonBoundAry;
			foreach my $exonID (keys %{$refExonRngByGeneHsh{$refTrnscptID}}) {
				push @refIsfmExonBoundAry, ${${$refExonRngByGeneHsh{$refTrnscptID}}{$exonID}}{"start"};
				push @refIsfmExonBoundAry, ${${$refExonRngByGeneHsh{$refTrnscptID}}{$exonID}}{"end"};
			}
			
			@refIsfmExonBoundAry = sort {$a<=>$b} @refIsfmExonBoundAry;
			my $geneStart = $refIsfmExonBoundAry[0];
			my $geneEnd = $refIsfmExonBoundAry[-1];

			my ($refDNASeq, $refFirst3Bases, $refLast3Bases, $refPepSeq, $refFullORF) = extractExonBoundSeqAndTranslate(\@refIsfmExonBoundAry, $cntgSeq, $strnd, \%selectedIDCodonTableHsh);

			#---generate GFF files
			my $refGenemRNALineStr = generatemRNAExonGFFLine($refTrnscptID, $refTrnscptID.".ref", \@refIsfmExonBoundAry, $cntg, $strnd, "y", $nameByGeneHsh{$refTrnscptID}, "completeORF=".$refFullORF, "y");
			my $refGeneOnlyLineStr = generatemRNAExonGFFLine($refTrnscptID, $refTrnscptID.".ref", \@refIsfmExonBoundAry, $cntg, $strnd, "y", $nameByGeneHsh{$refTrnscptID}, "completeORF=".$refFullORF, "n");

			print INDIVJREFGFF $refGenemRNALineStr."\n";
			print INDIVJBOTHGFF $refGenemRNALineStr."\n";
			print INDIVJALTGFF $refGeneOnlyLineStr."\n" if ($isofmNum > 0);

			#---goes through all alternative IndivJs
			foreach my $isofmID (keys %altIsfmJunctStrHsh) {

				my @altIsfmExonBoundAry;
				push @altIsfmExonBoundAry, $geneStart;
				push @altIsfmExonBoundAry, $geneEnd;
				foreach my $junctStr (keys %{$altIsfmJunctStrHsh{$isofmID}}) {
					my @junctStrSplt = split /:/, $junctStr;
					my $junctBoundStart = $junctStrSplt[1]-1;
					my $junctBoundEnd = $junctStrSplt[2]+1;
					push @altIsfmExonBoundAry, ($junctBoundStart, $junctBoundEnd);
				}
				@altIsfmExonBoundAry = sort {$a <=> $b} @altIsfmExonBoundAry;
				
				my ($altDNASeq, $altFirst3Bases, $altLast3Bases, $altPepSeq, $altFullORF) = extractExonBoundSeqAndTranslate(\@altIsfmExonBoundAry, $cntgSeq, $strnd, \%selectedIDCodonTableHsh);
				print FASTASEQ ">".$refTrnscptID."_".$isofmID."_".$altFullORF."\n";
				print FASTASEQ $altDNASeq."\n";
				my $NGSJunctStr = ${${$junctBasedAltIsofmInfoHsh{$refTrnscptID}}{$isofmID}}{"NGSJunctStr"};
				${$NGSJunctInfoHsh{$NGSJunctStr}}{"altIsofmFullORF"} = $altFullORF;

				my $altmRNAOnlyLineStr = generatemRNAExonGFFLine($refTrnscptID, $refTrnscptID.".".$isofmID, \@altIsfmExonBoundAry, $cntg, $strnd, "n", $nameByGeneHsh{$refTrnscptID}, "completeORF=".$altFullORF, "y");
				print INDIVJBOTHGFF $altmRNAOnlyLineStr."\n";
				print INDIVJALTGFF $altmRNAOnlyLineStr."\n";
			}
		}
	}
	close FASTASEQ;
	close INDIVJBOTHGFF;
	close INDIVJREFGFF;
	close INDIVJALTGFF;
	
	return \%altTrnscptInfoHsh, \%altJunctStrInfoHsh, \%refJunctInfoHsh, \%NGSJunctInfoHsh;

}
########################################################################## junctInfoToBEDLine
sub junctInfoToBEDLine {
	
	my $junctStr = $_[0];
	my $strnd =  $_[1];
	my $blk1Size =  $_[2];
	my $blk2Size = $_[3];
	my $blkName = $_[4];
	my $score = $_[5];
	
	my @junctStrSplt = split /\:/, $junctStr;
	my $cntg = $junctStrSplt[0];
	my $intronStart = $junctStrSplt[1];
	my $intronEnd = $junctStrSplt[2];
	
	my $bedStart = $intronStart - $blk1Size - 1;
	my $bedEnd = $intronEnd + $blk2Size;
	my $blk1Start = 0;
	my $blk2Start = $intronEnd - $intronStart + $blk2Size + 1;
	
	my $BEDLine = join "\t", ($cntg, $bedStart, $bedEnd, $blkName, $score, $strnd, $bedStart, $bedEnd, 0, 2, $blk1Size.",".$blk2Size.",", $blk1Start.",".$blk2Start.",");
	
	return $BEDLine;
	
}
########################################################################## generatemRNAExonGFFLine
sub generatemRNAExonGFFLine {

	my $geneParent = $_[0];
	my $mRNAID = $_[1];
	my @exonRng = @{$_[2]};
	my $cntg = $_[3];
	my $strnd = $_[4];
	my $outputGeneLine = $_[5];
	my $geneName = $_[6];
	my $mRNAOtherCommand = $_[7];
	my $outputmRNALine = $_[8];

	my @allOutputLineAry;
	
	#---assuming the mRNA and the gene are having the same range
	my $mRNAStart = $exonRng[0];
	my $mRNAEnd =  $exonRng[-1];
	my $geneLine = join "\t", ($cntg, "alternativeSplicer", "gene", $mRNAStart, $mRNAEnd, ".", $strnd, ".", "ID=".$geneParent.";"."Name=".$geneName.";");

	push @allOutputLineAry, $geneLine if ($outputGeneLine eq "y");
	if ($outputmRNALine eq "y") {
		my $mRNALine = join "\t", ($cntg, "alternativeSplicer", "mRNA", $mRNAStart, $mRNAEnd, ".", $strnd, ".", "ID=".$mRNAID.";"."Name=".$mRNAID.";"."Parent=".$geneParent.";".$mRNAOtherCommand);
		push @allOutputLineAry, $mRNALine;
		my $exonNum = 0;
		for (my $i = 0; $i < $#exonRng; $i += 2) {#---$i = 0, 2, 4, 6 is there're are 8 elements
			my $exonStart = $exonRng[$i];
			my $exonEnd = $exonRng[$i+1];
			$exonNum++;
			my $exonLine = join "\t", ($cntg, "alternativeSplicer", "exon", $exonStart, $exonEnd, ".", $strnd, ".", "ID=exon_".$mRNAID."-".$exonNum.";"."Name=exon;"."Parent=".$mRNAID.";");
			push @allOutputLineAry, $exonLine;
		}
	}
	
	my $allOutputLineStr = join "\n", @allOutputLineAry;
	
	return $allOutputLineStr;
}
########################################################################## readMultiFasta
sub readMultiFasta {

	my $refFastaPath = $_[0];
	my ($seq, $seqName, %fastaHsh);
	my $i = 0;
	print "Reading $refFastaPath into a hash.\n";
	open (INFILE, $refFastaPath);
	chomp (my $curntLine = <INFILE>); #get the first line
	while (my $nextLine = <INFILE>) {
		chomp $nextLine;
		
		#---Only two types of line in current line, the header or seq
		if ($curntLine =~ m/^>/) {#-- header line
			my @theLineSplt = split (/\|/, $curntLine);
			$seqName = $theLineSplt[0]; #---get the first tag
			$seqName =~ s/ //g; #---remove space
			$seqName =~ s/>//g; #---remove space
		} else {#--seq line
			$seq = $seq.$curntLine;
		}
		
		#---check if next line has a > or that's the end of file
		if ($nextLine =~ m/^>/) {
			$fastaHsh{$seqName} = $seq;
			$seq = "";
		} elsif (eof(INFILE)) {#---this is the last line
			$seq = $seq.$nextLine;
			$fastaHsh{$seqName} = $seq;
		}
		
		#---next line becomes current line
		$curntLine = $nextLine;
	}
	close INFILE;
	
	return (\%fastaHsh);
}
########################################################################## extractExonBoundSeqAndTranslate
sub extractExonBoundSeqAndTranslate {

	my @tmpExonBoundAry = @{$_[0]};
	my $cntgSeq = $_[1];
	my $strnd = $_[2];
	my %selectedIDCodonTableHsh = %{$_[3]};

	#---check if the bounds are in pair
	die "tmpExonBoundAry must contains even number elements\n" if (@tmpExonBoundAry/2 != (int  @tmpExonBoundAry/2)); 
	
	my $DNASeqToReturn = "";
	
	for (my $i = 0; $i < $#tmpExonBoundAry; $i += 2) {#---$i = 0, 2, 4, 6 is there're are 8 elements
		my $exonStart = $tmpExonBoundAry[$i];
		my $exonEnd = $tmpExonBoundAry[$i+1];
		my $length = $exonEnd-$exonStart+1;
		my $offset = $exonStart-1;
		$DNASeqToReturn .= substr ($cntgSeq, $offset, $length);
	}
	
	if ($strnd eq "-") {
		$DNASeqToReturn = reverse $DNASeqToReturn;
		$DNASeqToReturn =~ tr/ACGTacgt/TGCAtgca/;
	}
	
	my @codonAry = $DNASeqToReturn =~ /\w{3}/g;
	
	my $pepSeqToReturn = "";
	foreach my $codon (@codonAry) {
		my $aminoAcid = "X"; 
		$aminoAcid = $selectedIDCodonTableHsh{$codon} if (exists $selectedIDCodonTableHsh{$codon});
		$pepSeqToReturn .= $aminoAcid;
	}
	
	my @DNASeqToReturnSplt = split //, $DNASeqToReturn;
	
	my $first3Bases = $codonAry[0];
	my $last3Bases = $codonAry[-1];
	#print "$first3Bases $last3Bases $strnd                     \n";
	
	chop $pepSeqToReturn; #---dont not return the stop codon
	
	my $fullORF = "y";
	$fullORF = "n" if (($pepSeqToReturn =~ m/X/) or ($pepSeqToReturn =~ m/\*/));
	
	return ($DNASeqToReturn, $first3Bases, $last3Bases, $pepSeqToReturn, $fullORF);
}
########################################################################## getCodonTable
sub getCodonTable {
	
	my $tablePath = $_[0];
	my $selectedID = $_[1]; #---use all to get all ID, or a number of get a specific one
	my $ID = "none";
	my $name = "none";
	my %codonTableStrByIDHsh;
	open (INFILE, $tablePath);
	while (my $theLine = <INFILE>) {
		chomp $theLine;
		$name = "" if ($theLine =~ m/^ \{/);
		$name .= $1.";" if ($theLine =~ m/^ +name +\"(.+)\"/);
		if ($theLine =~ m/^ +id +([0-9]+) +,$/) {
			$ID = $1;
			${$codonTableStrByIDHsh{$ID}}{"name"} = $name;
		}
		${$codonTableStrByIDHsh{$ID}}{"aminoAcid"} = $1 if ($theLine =~ m/^ +ncbieaa +\"(.+)\"/);
		${$codonTableStrByIDHsh{$ID}}{"altStart"} = $1 if ($theLine =~ m/^ +sncbieaa +\"(.+)\"/);
		${$codonTableStrByIDHsh{$ID}}{"base1"} = $1 if ($theLine =~ m/^ +-- +Base1 +(.+)$/);
		${$codonTableStrByIDHsh{$ID}}{"base2"} = $1 if ($theLine =~ m/^ +-- +Base2 +(.+)$/);
		${$codonTableStrByIDHsh{$ID}}{"base3"} = $1 if ($theLine =~ m/^ +-- +Base3 +(.+)$/);
	}
	close INFILE;
	
	my (%allIDCodonTableHsh, %allIDStartCodonHsh, %selectedIDCodonTableHsh, %selectedIDStartCodonHsh, @tmpStopCodonAry, @tmpStartCodonAry, @tmpAltStartAty);
	
	foreach my $ID (sort {$a <=>$b} keys %codonTableStrByIDHsh) {
		my @aminoAcidAry = split '', ${$codonTableStrByIDHsh{$ID}}{"aminoAcid"}; 		
		my @altStartAry = split '', ${$codonTableStrByIDHsh{$ID}}{"altStart"}; 		
		my @base1Ary = split '', ${$codonTableStrByIDHsh{$ID}}{"base1"}; 		
		my @base2Ary = split '', ${$codonTableStrByIDHsh{$ID}}{"base2"};
		my @base3Ary = split '', ${$codonTableStrByIDHsh{$ID}}{"base3"};
		
		for (my $i = 0; $i < @base1Ary; $i++ ) {
			my $theCodon = $base1Ary[$i].$base2Ary[$i].$base3Ary[$i];
			my $theAminoAcid = $aminoAcidAry[$i];
			my $altStart = $altStartAry[$i];
			
			$altStart = "-" if (($ATGOnly eq "y") and ($theCodon ne "ATG"));
			
			${$allIDCodonTableHsh{$ID}}{$theCodon} = $theAminoAcid;
			${$allIDStartCodonHsh{$ID}}{$theCodon} = $altStart;
						
			if ($ID eq $selectedID) {
				$selectedIDCodonTableHsh{$theCodon} = $theAminoAcid;
				$selectedIDStartCodonHsh{$theCodon} = $altStart;
				push @tmpAltStartAty, $altStart;
				push @tmpStopCodonAry, $theCodon if ($theAminoAcid eq "*"); 
				push @tmpStartCodonAry, $theCodon if ($altStart eq "M"); 
			}
		}
	}
	
	if ($selectedID eq "all") {
		return (\%allIDCodonTableHsh, \%allIDStartCodonHsh);
	} else {
		my $tmp = keys %selectedIDStartCodonHsh;
		my $allStartCodon = join ";", @tmpStartCodonAry;
		my $allStopCodon = join ";", @tmpStopCodonAry;
		die "ID $selectedID is not found in the supplied codon table" if ($tmp == 0);
		print "You have chosen ".${$codonTableStrByIDHsh{$selectedID}}{"name"}." genetic code\n\n";
		print "Start codon = $allStartCodon\n";
		print "Stop codon = $allStopCodon\n\n";

		${$codonTableStrByIDHsh{$selectedID}}{"altStart"} = join "", @tmpAltStartAty;
		print "The Codon Table with ATGOnly=$ATGOnly:\n";
		print ${$codonTableStrByIDHsh{$selectedID}}{"aminoAcid"}."\n";
		print ${$codonTableStrByIDHsh{$selectedID}}{"altStart"}."\n"; 
		print ${$codonTableStrByIDHsh{$selectedID}}{"base1"}."\n";
		print ${$codonTableStrByIDHsh{$selectedID}}{"base2"}."\n";
		print ${$codonTableStrByIDHsh{$selectedID}}{"base3"}."\n\n";
		return (\%selectedIDCodonTableHsh, \%selectedIDStartCodonHsh);
	}
}
########################################################################## getPileUpInfo
sub getUnqAndLoComInfo {
	
	my %cntgSeqHsh = %{$_[0]};
	my $unqPileupPath = $_[1];
	my $lowComRegionPath = $_[2];
	
	#---create an emtpy 
	my %unqLoComByCntgPosHsh;

	if (($lowComRegionPath ne "no") or ($unqPileupPath ne "no")) {
		printProgressScale("Creating a hash to contain the pileup info", 50);
		my $totalCntg = keys %cntgSeqHsh;
		my $procCntg = 0;
		foreach my $cntg (sort {$a cmp $b} keys %cntgSeqHsh) {
			$procCntg++;
			updateProgressBar($cntg, $procCntg, $totalCntg, 50, 10);
			my $cntgLength = length $cntgSeqHsh{$cntg};
			foreach my $pos (1..$cntgLength) {
				push @{$unqLoComByCntgPosHsh{$cntg}}, "0"."\t"."0"."\t"."0";
			}
		}
	}

	if ($unqPileupPath ne "no") {
		open (INFILE, "$unqPileupPath");
		my ($fileToCheckSizeTotalLineNum, $intervalSize) = checkFileSizeAndDefineIntervalSize($unqPileupPath, 10000);
		my $intervalStart = time(); 
		my $lineProc = my $progCount = 0;
		while (my $theLine = <INFILE>) {
			$lineProc++;$progCount++;
			($progCount, $intervalStart) = reportProgress($progCount, $lineProc, $intervalSize, $fileToCheckSizeTotalLineNum, $intervalStart) if ($progCount == $intervalSize);
			my @theLineSplt = split /\t/, $theLine;
			my $cntg = $theLineSplt[0];
			my $pos = $theLineSplt[1];
			my $futr = $theLineSplt[2];
			my $plusUnq = $theLineSplt[5];
			my $minusUnq = $theLineSplt[6];
			my $index = $pos-1;
			my @unqLoComSplt = split /\t/, ${$unqLoComByCntgPosHsh{$cntg}}[$index];
			$unqLoComSplt[0] = $plusUnq;
			$unqLoComSplt[1] = $minusUnq;
			${$unqLoComByCntgPosHsh{$cntg}}[$index] = join "\t", @unqLoComSplt;
		}
		close INFILE;
	} else {
		
		print "\nSkip reading pileup file\n";
		
	}
	
	if ($lowComRegionPath ne "no") {
		print "Reading low complexity region\n";
		
		open (LOWCOMFILE, "$lowComRegionPath");
		my ($fileToCheckSizeTotalLineNum, $intervalSize) = checkFileSizeAndDefineIntervalSize($lowComRegionPath, 10000);
		my $intervalStart = time(); 
		my $lineProc = my $progCount = 0;
		my $cntg = "null";
		while (my $theLine = <LOWCOMFILE>) {
			chomp $theLine;
			$lineProc++;$progCount++;
			($progCount, $intervalStart) = reportProgress($progCount, $lineProc, $intervalSize, $fileToCheckSizeTotalLineNum, $intervalStart) if ($progCount == $intervalSize);
			$cntg = substr ($theLine, index ($theLine, ">")+1) if ($theLine =~ m/^>/);
			if ($theLine =~ m/^\d+? - \d+?$/) {
				my @theLineSplt = split / - /, $theLine;
				my $loCmStart = $theLineSplt[0];
				my $loCmEnd = $theLineSplt[1];
				for my $pos ($loCmStart..$loCmEnd) {
					my $index = $pos - 1;
					my @unqLoComSplt = split /\t/, ${$unqLoComByCntgPosHsh{$cntg}}[$index];
					$unqLoComSplt[2] = 1;
					${$unqLoComByCntgPosHsh{$cntg}}[$index] = join "\t", @unqLoComSplt;
				}
			}
		}
		close LOWCOMFILE;
	} else {
		
		print "\nSkip reading lowComRegion file\n";
		
	}

	return \%unqLoComByCntgPosHsh;
}
########################################################################## checkFileSizeAndDefineIntervalSize
sub checkFileSizeAndDefineIntervalSize {
    
    my $fileToCheckPath = $_[0];
    my $linesToSample = $_[1];
	
	#---make sure $linesToSample is a non-zero number, if not set to 10000
	my $linesToSampleInt = int $linesToSample;
	$linesToSampleInt = 100000 if (($linesToSampleInt != $linesToSample) or ($linesToSampleInt == 0));

    #---get the filename from the path
    my @fileToCheckPathSplt = split /\//, $fileToCheckPath;
    my $fileToCheckName = $fileToCheckPathSplt[-1];
    
    print "Estimating the number of lines in $fileToCheckName.\n";
	
	#---estimate the number of lines in the file
	open (INFILE, $fileToCheckPath) || die "Can't open $fileToCheckPath.\n";
	my $tmpFilePath = $fileToCheckPath."_tmp.txt";
	system "tail -$linesToSampleInt $fileToCheckPath >$tmpFilePath";
	my $fileToCheckSize = -s "$fileToCheckPath";
   	my $tmpFileSize = -s "$tmpFilePath";
	system "rm $tmpFilePath";
   	my $fileToCheckSizeTotalLineNum = int (($fileToCheckSize/$tmpFileSize)*$linesToSampleInt);
   	print "Estimated to have ".$fileToCheckSizeTotalLineNum." lines in $fileToCheckName.\n";

	my $intervalSize = int ($fileToCheckSizeTotalLineNum/100); #---define as 
	$intervalSize = 1000000 if  ($intervalSize > 1000000);
	
	return ($fileToCheckSizeTotalLineNum, $intervalSize);

}
########################################################################## reportProgress
sub reportProgress {

	my $progCount = $_[0];
	my $lineProc = $_[1];
	my $intervalSize = $_[2];
	my $fileTotalLineNum = $_[3];
	my $intervalStart = $_[4];

	$progCount=0;
	my $intervalEnd = time();
	my $timeElapsed = $intervalEnd - $intervalStart;
	$timeElapsed = sprintf ("%.2f", $timeElapsed);
	my $estimatedEnd = (($fileTotalLineNum - $lineProc)*$timeElapsed)/$intervalSize;
	$estimatedEnd = sprintf ("%.2f", $estimatedEnd/60);
	print "$lineProc lines processed. Last $intervalSize lines:".$timeElapsed." sec. Estimated end: ".$estimatedEnd." mins.\r";
	$intervalStart = time();
	
	return ($progCount, $intervalStart);
		
}
########################################################################## getJunctUnqAndSeq
sub getJunctUnqAndSeq {
	
	my %NGSJunctCntgHsh = %{$_[0]};
	my %unqLoComByCntgPosHsh = %{$_[1]};
	my %NGSJunctStrndHsh = %{$_[2]};
	my %fastaHsh = %{$_[3]};
	my %refGeneIDByJunctStrHsh = %{$_[4]};
	my %refStrndHsh = %{$_[5]};
	my %SSOvrlpRefJNGSJHsh = %{$_[6]};
	my %nonCanNGSJunctInfoHsh = %{$_[7]};
	my $exonSeqLength = $_[8];
	
	#---covert to hash keys by Cntg
	my (%junctStrByCntgHsh, %allJunctStrndHsh);
	
	my (%NGSJunctUnqHsh, %NGSJunctSeqHsh, %refJunctInfoHsh);
	
	foreach my $junctStr (keys %nonCanNGSJunctInfoHsh) {
		my @junctStrSplt = split /:/, $junctStr;
		my $cntg = $junctStrSplt[0];
		${$junctStrByCntgHsh{$cntg}}{$junctStr}++;
		my $strnd = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"strnd"};
		$strnd = "+" if ($strnd eq "*"); #---if no strand, treat it as +
		$allJunctStrndHsh{$junctStr} = $strnd;
		${$junctStrByCntgHsh{$cntg}}{$junctStr}++;
	}
	
	foreach my $junctStr (keys %NGSJunctCntgHsh) {
		my @junctStrSplt = split /:/, $junctStr;
		my $cntg = $junctStrSplt[0];
		${$junctStrByCntgHsh{$cntg}}{$junctStr}++;
		$allJunctStrndHsh{$junctStr} = $NGSJunctStrndHsh{$junctStr};
	}

	foreach my $junctStr (keys %refGeneIDByJunctStrHsh) {
		my @junctStrSplt = split /:/, $junctStr;
		my $cntg = $junctStrSplt[0];
		my $geneID = $refGeneIDByJunctStrHsh{$junctStr};
		${$refJunctInfoHsh{$junctStr}}{"geneID"} = $geneID;
		${$refJunctInfoHsh{$junctStr}}{"strnd"} = $refStrndHsh{$geneID};
		${$refJunctInfoHsh{$junctStr}}{"cntg"} = $cntg;
		${$junctStrByCntgHsh{$cntg}}{$junctStr}++;
		$allJunctStrndHsh{$junctStr} = $refStrndHsh{$geneID};
		
		${$refJunctInfoHsh{$junctStr}}{"cntg"} = $cntg;
		my ($junctCntg, $junctStart, $junctStop) = split /:/, $junctStr;
		my $intronSize = $junctStop - $junctStart + 1;
		my $intronFrame = '';
				
		if ($intronSize/3 == int ($intronSize/3)) {
			$intronFrame = '3n';
		} elsif (($intronSize+1)/3 == int (($intronSize+1)/3)) {
			$intronFrame = '3n+1';
		
		} elsif (($intronSize+2)/3 == int (($intronSize+2)/3)) {
			$intronFrame = '3n+2';
		
		} else {
			die "Impossible intronFrame\n";
		}
				
		${$refJunctInfoHsh{$junctStr}}{"intronSize"} = $intronSize;
		${$refJunctInfoHsh{$junctStr}}{"intronFrame"} = $intronFrame;
	}
	
	my $totalJunctNum = keys %allJunctStrndHsh;
	my $procJunctNum = 0;
	
	printProgressScale("\nGetting unq for NGS junctions", 50);
	
	#----read the precalculated unq file
	my (%tmpJunctUnqHsh, %tmpJunctLoCmHsh);
	if (($unqPileupPath eq "no") or ($lowComRegionPath eq "no")) {
		open PRECALUNQFILE, "$precalUnqLoCmPath";
		while (my $theLine = <PRECALUNQFILE>) {
			chomp $theLine;
			my @theLineSplt = split /\t/, $theLine;
			$tmpJunctUnqHsh{$theLineSplt[0]} = $theLineSplt[1];
			$tmpJunctLoCmHsh{$theLineSplt[0]} = $theLineSplt[2];
		}
		close PRECALUNQFILE;
	}
	
	open OUTUNQFILE, ">$outDir/preCalInfo/precalculated.unq.txt" if (($unqPileupPath ne "no") or ($lowComRegionPath ne "no"));
	
	#---go through all cntg and all junctStr
	foreach my $cntg (keys %junctStrByCntgHsh) {
		my @tmpCntgUnqAry;
		if (($unqPileupPath ne "no") or ($lowComRegionPath ne "no")) {
			@tmpCntgUnqAry = @{$unqLoComByCntgPosHsh{$cntg}}; #---get the cntg unq
		}
		my $cntgSeq = $fastaHsh{$cntg};
		my $endIndex = $#tmpCntgUnqAry;
		
		foreach my $junctStr (keys %{$junctStrByCntgHsh{$cntg}}) {
			$procJunctNum++;
			updateProgressBar($junctStr, $procJunctNum, $totalJunctNum, 50, 10);

			#---get the seq
			my @junctStrSplt = split /:/, $junctStr;
			my $junctStart = $junctStrSplt[1]-1;
			my $junctEnd = $junctStrSplt[2]-1;
			my $upStrmBound = $junctStart - $boundWidth;
			my $downStrmBound = $junctEnd + $boundWidth;
			$upStrmBound = 1 if ($upStrmBound < 1);
			$downStrmBound = $endIndex if ($downStrmBound > $endIndex);
			my @boundRngAry = ($upStrmBound, $junctStart, $junctEnd, $downStrmBound);
			my $strnd = $allJunctStrndHsh{$junctStr};

			my $leftFlankSeq = substr $cntgSeq, $junctStart-$boundWidth, ($boundWidth*2)+2;
			my $rightFlankSeq = substr $cntgSeq, $junctEnd-$boundWidth-1, ($boundWidth*2)+2;
			my $left2ntJnctnSeq = substr $cntgSeq, $junctStart, 2;
			my $right2ntJnctnSeq = substr $cntgSeq, $junctEnd-1, 2;
			my $leftExonSeq = substr $cntgSeq, (($junctStart-1)-$exonSeqLength), $exonSeqLength;
			my $rightExonSeq = substr $cntgSeq, $junctEnd+1, $exonSeqLength;
			
			
			my $completeSeq = "complete";
			if ((length $leftFlankSeq != ($boundWidth*2)+2) or (length $rightFlankSeq != ($boundWidth*2)+2)) {
				print "WARNING: $junctStr is too close to contig end.\r";
				$completeSeq = "cntgEnd";
			}

			if (($leftFlankSeq =~ m/N/gi) or ($rightFlankSeq =~ m/N/gi)) {
				print "WARNING: $junctStr is too close to N regions.\r";
				$completeSeq = "NReg";
			}
			
			if ($strnd eq "+") {

			} elsif ($strnd eq "-") {
				($left2ntJnctnSeq, $right2ntJnctnSeq) = ($right2ntJnctnSeq, $left2ntJnctnSeq);
				($leftFlankSeq, $rightFlankSeq) = ($rightFlankSeq, $leftFlankSeq);
				($leftExonSeq, $rightExonSeq) = ($rightExonSeq, $leftExonSeq);
				
				$leftFlankSeq = reverse $leftFlankSeq;
				$rightFlankSeq = reverse $rightFlankSeq;
				$leftFlankSeq =~ tr/ACGTacgt/TGCAtgca/;
				$rightFlankSeq =~ tr/ACGTacgt/TGCAtgca/;

				$left2ntJnctnSeq = reverse $left2ntJnctnSeq;
				$right2ntJnctnSeq = reverse $right2ntJnctnSeq;
				$left2ntJnctnSeq =~ tr/ACGTacgt/TGCAtgca/;
				$right2ntJnctnSeq =~ tr/ACGTacgt/TGCAtgca/;

				$leftExonSeq = reverse $leftExonSeq;
				$rightExonSeq = reverse $rightExonSeq;
				$leftExonSeq =~ tr/ACGTacgt/TGCAtgca/;
				$rightExonSeq =~ tr/ACGTacgt/TGCAtgca/;

			} else {
				print "WARNING: undefined strand info.\n";
			}
			
			$leftFlankSeq =~ tr/Tt/Uu/;
			$rightFlankSeq =~ tr/Tt/Uu/;
			$left2ntJnctnSeq =~ tr/Tt/Uu/;
			$right2ntJnctnSeq =~ tr/Tt/Uu/;
			$leftExonSeq =~ tr/Tt/Uu/;
			$rightExonSeq =~ tr/Tt/Uu/;
			my $leftExonPolyBasePrptn;
			my $rightExonPolyBasePrptn;
			my @baseToCountAry = ("A", "U", "G", "C", "N");
			($dummy, $leftExonPolyBasePrptn, $dummy, $dummy, $dummy) = countBaseContent($leftExonSeq, \@baseToCountAry);
			($dummy, $rightExonPolyBasePrptn, $dummy, $dummy, $dummy) = countBaseContent($rightExonSeq, \@baseToCountAry);

			#print $leftExonSeq."\t".$leftExonPolyBasePrptn."\n";
			#print $rightExonSeq."\t".$rightExonPolyBasePrptn."\n";
			
			#---get the unq
			my $unq;
			my $loCmPrptn;
			if (($unqPileupPath ne "no") or ($lowComRegionPath ne "no")) {
				my $unqSum = my $nonZeroUnqPosNum = my $loCmSum = my $loCmPos = 0;
				for (my $i = 0; $i < $#boundRngAry; $i += 2) {#---$i = 0, 2, 4, 6 is there're are 8 elements
					my $startPos = $boundRngAry[$i]-1;
					my $endPos = $boundRngAry[$i+1]-1;
					
					foreach my $unqStr (@tmpCntgUnqAry[$startPos..$endPos]) {
						my @unqStrSplt = split /\t/, $unqStr;
						my $plusUnq = $unqStrSplt[0];
						my $minusUnq = $unqStrSplt[1];
						my $loCmReg = $unqStrSplt[2];

						my $strndUnq;
						
						if ($unqPileupPath ne "no") {
							if ($strnd eq "+") {
								$strndUnq = $plusUnq;
							} elsif ($strnd eq "-") {
								$strndUnq = $minusUnq;
							} else {
								print "WARNING: undefined strand info.\n";
							}
							
							if ($strndUnq > 0) {
								$unqSum += $strndUnq;
								$nonZeroUnqPosNum++;
							}
						}
						
						if ($lowComRegionPath ne "no") {
							$loCmSum += $loCmReg;
							$loCmPos++;
						}
					}
				}

				if ($nonZeroUnqPosNum < ($boundWidth/2)) {#---defined as undefined is less than half of $boundWidth pos counted, on cntg edge
					$unq = -1;
				} else {
					$unq = sprintf '%.06f', $unqSum/$nonZeroUnqPosNum;
				}

				$loCmPrptn = sprintf '%.06f', $loCmSum/$loCmPos;
			}

			if ($unqPileupPath eq "no") {
				if (exists $tmpJunctUnqHsh{$junctStr}) {
					$unq = $tmpJunctUnqHsh{$junctStr};
				} else {
					print "warning: unq of $junctStr is not defined. set to -1\n";
					$unq = -1;
				}
			}

			if ($lowComRegionPath eq "no") {
				if (exists $tmpJunctLoCmHsh{$junctStr}) {
					$loCmPrptn = $tmpJunctLoCmHsh{$junctStr};
				} else {
					print "warning: loCmPrptn of $junctStr is not defined. set to -1\n";
					$loCmPrptn = -1;
				}
			}
			
			print OUTUNQFILE $junctStr."\t".$unq."\t".$loCmPrptn."\n"  if (($unqPileupPath ne "no") or ($lowComRegionPath ne "no"));;
			
			#---store the info, if it is an NGS Junct
			if (exists $NGSJunctCntgHsh{$junctStr}) {
				$NGSJunctUnqHsh{$junctStr} = $unq;
				${$NGSJunctSeqHsh{$junctStr}}{"loCmPrptn"} = $loCmPrptn;
				${$NGSJunctSeqHsh{$junctStr}}{"unq"} = $unq;
				${$NGSJunctSeqHsh{$junctStr}}{"leftFlankSeq"} = $leftFlankSeq;
				${$NGSJunctSeqHsh{$junctStr}}{"rightFlankSeq"} = $rightFlankSeq;
				${$NGSJunctSeqHsh{$junctStr}}{"upSite"} = $left2ntJnctnSeq;
				${$NGSJunctSeqHsh{$junctStr}}{"downSite"} = $right2ntJnctnSeq;
				${$NGSJunctSeqHsh{$junctStr}}{"completeSeq"} = $completeSeq;
				${$NGSJunctSeqHsh{$junctStr}}{"leftExonSeq"} = $leftExonSeq;
				${$NGSJunctSeqHsh{$junctStr}}{"rightExonSeq"} = $rightExonSeq;
				${$NGSJunctSeqHsh{$junctStr}}{"leftExonPolyBasePrptn"} = $leftExonPolyBasePrptn;
				${$NGSJunctSeqHsh{$junctStr}}{"rightExonPolyBasePrptn"} = $rightExonPolyBasePrptn;
			}
			
			#---store the info, if it is an nonCan NGS Junct
			if (exists $nonCanNGSJunctInfoHsh{$junctStr}) {
				${$nonCanNGSJunctInfoHsh{$junctStr}}{"loCmPrptn"} = $loCmPrptn;
				${$nonCanNGSJunctInfoHsh{$junctStr}}{"unq"} = $unq;
				${$nonCanNGSJunctInfoHsh{$junctStr}}{"leftFlankSeq"} = $leftFlankSeq;
				${$nonCanNGSJunctInfoHsh{$junctStr}}{"rightFlankSeq"} = $rightFlankSeq;
				${$nonCanNGSJunctInfoHsh{$junctStr}}{"upSite"} = $left2ntJnctnSeq;
				${$nonCanNGSJunctInfoHsh{$junctStr}}{"downSite"} = $right2ntJnctnSeq;
				${$nonCanNGSJunctInfoHsh{$junctStr}}{"completeSeq"} = $completeSeq;
				${$nonCanNGSJunctInfoHsh{$junctStr}}{"leftExonSeq"} = $leftExonSeq;
				${$nonCanNGSJunctInfoHsh{$junctStr}}{"rightExonSeq"} = $rightExonSeq;
				${$nonCanNGSJunctInfoHsh{$junctStr}}{"leftExonPolyBasePrptn"} = $leftExonPolyBasePrptn;
				${$nonCanNGSJunctInfoHsh{$junctStr}}{"rightExonPolyBasePrptn"} = $rightExonPolyBasePrptn;
			}
			
			#---store the info, if it is an ref Junct
			if (exists $refGeneIDByJunctStrHsh{$junctStr}) {

				${$refJunctInfoHsh{$junctStr}}{"loCmPrptn"} = $loCmPrptn;
				${$refJunctInfoHsh{$junctStr}}{"NGSHit"} = "missed";
				${$refJunctInfoHsh{$junctStr}}{"unq"}= $unq;
				${$refJunctInfoHsh{$junctStr}}{"strnd"}= $strnd;
				${$refJunctInfoHsh{$junctStr}}{"leftFlankSeq"} = $leftFlankSeq;
				${$refJunctInfoHsh{$junctStr}}{"rightFlankSeq"} = $rightFlankSeq;
				${$refJunctInfoHsh{$junctStr}}{"upSite"} = $left2ntJnctnSeq;
				${$refJunctInfoHsh{$junctStr}}{"downSite"} = $right2ntJnctnSeq;
				${$refJunctInfoHsh{$junctStr}}{"leftExonSeq"} = $leftExonSeq;
				${$refJunctInfoHsh{$junctStr}}{"rightExonSeq"} = $rightExonSeq;
				${$refJunctInfoHsh{$junctStr}}{"leftExonPolyBasePrptn"} = $leftExonPolyBasePrptn;
				${$refJunctInfoHsh{$junctStr}}{"rightExonPolyBasePrptn"} = $rightExonPolyBasePrptn;
				
				if (exists $SSOvrlpRefJNGSJHsh{$junctStr}) {
					${$refJunctInfoHsh{$junctStr}}{"NGSHit"} = "partialOverlap";
					foreach my $NGSJunctStr (keys %{$SSOvrlpRefJNGSJHsh{$junctStr}}) {
						${$refJunctInfoHsh{$junctStr}}{"NGSHit"} = "exactOverlap" if ($junctStr eq $NGSJunctStr);
					}
				}
				${$refJunctInfoHsh{$junctStr}}{"completeSeq"} = $completeSeq;
			}
		}
	}
	
	print "\n";
	
	return (\%NGSJunctUnqHsh, \%NGSJunctSeqHsh, \%refJunctInfoHsh, \%nonCanNGSJunctInfoHsh);
}
########################################################################## sortAndPrintAllJunctInfo
sub sortAndPrintAllJunctInfo {
	
	my %NGSJunctInfoHsh = %{$_[0]};
	my %refJunctInfoHsh = %{$_[1]};
	my %refJunctCommentHsh = %{$_[2]};
	my %jClusterInfoHsh = %{$_[3]};
	
	my $trimExonEnd = 40;
	my $trimIntronEnd = 40;
	
	my %intronFrameHsh; #---to contain the intronFrameInfo
	
	open ALLNGSJLOG, ">$outDir/junctInfo/all_NGS_junct_log.txt";
	print ALLNGSJLOG join "", ((join "\t", ("junctStr", "readNum", "score", "strnd", "unq", "upSite", "downSite", "cluster", "majorCluster", "clstJunctNum", "prmntReadNumRatio", "majorJunct", "prmntInCluster", "onPrmntIsofm", "superJunct", "hitRef", "completeSeq", "leftFlankSeq", "rightFlankSeq", "splicingRatio", "BEDCoverRng", "ovrlpCtrgy", "ovrlpSense", "ovrlpAntisense", "superCtgry", "senseSplicingEfficiency", "antisenseSplicingEfficiency", "prmntDiff5", "prmntDiff3", "prmntDiffBoth", "splicingEffType", "geneCovNt", "refOrAlt", "isofmID", "altIsofmFullORF", "altSplcType", "onRefGene", "refJunctAvgReadNum", "NGSRefJunctReadNumRatio", "readNumCovPerNtRatio", "altSiteType", "exonSkipType", "refShift5", "refShift3", "skippedExonNum", "ovrlpRefJunctRdRatio", "avgOvrlpRefJReadNum", "consensusValue", "intronSize", "intronFrame")), "\n");

	my (%allSeqPathByJunctTypeHsh, %junctTypeCountHsh, %junctTypeScoreHsh, %splicingSiteDiffAllHsh);
	$allSeqPathByJunctTypeHsh{"NGSMajorClstPrmnt"} = "$outDir/flankSeq/NGSMajorClstPrmnt.fasta";
	$allSeqPathByJunctTypeHsh{"NGSMinorNonClstPrmnt"} = "$outDir/flankSeq/NGSMinorNonClstPrmnt.fasta";
	$allSeqPathByJunctTypeHsh{"NGSMajor"} = "$outDir/flankSeq/NGSMajor.fasta";
	$allSeqPathByJunctTypeHsh{"NGSMinor"} = "$outDir/flankSeq/NGSMinor.fasta";
	$allSeqPathByJunctTypeHsh{"NGSClstPrmnt"} = "$outDir/flankSeq/NGSClstPrmnt.fasta";
	$allSeqPathByJunctTypeHsh{"NGSNonClstPrmnt"} = "$outDir/flankSeq/NGSNonClstPrmnt.fasta";
	$allSeqPathByJunctTypeHsh{"NGSInfMajorClstPrmnt"} = "$outDir/flankSeq/NGSInfMajorClstPrmnt.fasta";
	$allSeqPathByJunctTypeHsh{"NGSInfMajorNonClstPrmnt"} = "$outDir/flankSeq/NGSInfMajorNonClstPrmnt.fasta";
	$allSeqPathByJunctTypeHsh{"NGSInfMajorNonClstPrmntLimA"} = "$outDir/flankSeq/NGSInfMajorNonClstPrmntLimA.fasta";#---limited to certain ratio to Prmnt J
	$allSeqPathByJunctTypeHsh{"NGSInfMajorNonClstPrmntLimB"} = "$outDir/flankSeq/NGSInfMajorNonClstPrmntLimB.fasta";#---limited to certain ratio to Prmnt J
	$allSeqPathByJunctTypeHsh{"NGSInfMajorNonClstPrmntLimC"} = "$outDir/flankSeq/NGSInfMajorNonClstPrmntLimC.fasta";#---limited to certain ratio to Prmnt J
	$allSeqPathByJunctTypeHsh{"NGSInfMajorNonClstPrmntLimD"} = "$outDir/flankSeq/NGSInfMajorNonClstPrmntLimD.fasta";#---limited to certain ratio to Prmnt J
	$allSeqPathByJunctTypeHsh{"NGSInfMajorNonClstPrmntLimE"} = "$outDir/flankSeq/NGSInfMajorNonClstPrmntLimE.fasta";#---limited to certain ratio to Prmnt J

	$allSeqPathByJunctTypeHsh{"NGSSplicingRatioLimA"} = "$outDir/flankSeq/NGSSplicingRatioLimA.fasta";
	$allSeqPathByJunctTypeHsh{"NGSSplicingRatioLimB"} = "$outDir/flankSeq/NGSSplicingRatioLimB.fasta";
	$allSeqPathByJunctTypeHsh{"NGSSplicingRatioLimC"} = "$outDir/flankSeq/NGSSplicingRatioLimC.fasta";
	$allSeqPathByJunctTypeHsh{"NGSSplicingRatioLimD"} = "$outDir/flankSeq/NGSSplicingRatioLimD.fasta";
	$allSeqPathByJunctTypeHsh{"NGSSplicingRatioLimE"} = "$outDir/flankSeq/NGSSplicingRatioLimE.fasta";
	
	$allSeqPathByJunctTypeHsh{"NGSScoreLimA"} = "$outDir/flankSeq/NGSScoreLimA.fasta";
	$allSeqPathByJunctTypeHsh{"NGSScoreLimB"} = "$outDir/flankSeq/NGSScoreLimB.fasta";
	$allSeqPathByJunctTypeHsh{"NGSScoreLimC"} = "$outDir/flankSeq/NGSScoreLimC.fasta";
	$allSeqPathByJunctTypeHsh{"NGSScoreLimD"} = "$outDir/flankSeq/NGSScoreLimD.fasta";
	$allSeqPathByJunctTypeHsh{"NGSScoreLimE"} = "$outDir/flankSeq/NGSScoreLimE.fasta";
	
	$allSeqPathByJunctTypeHsh{"NGSRefHit"} = "$outDir/flankSeq/NGSRefHit.fasta";
	$allSeqPathByJunctTypeHsh{"NGSNonRefHitGreater1000X"} = "$outDir/flankSeq/NGSNonRefHitGreater1000X.fasta";
	$allSeqPathByJunctTypeHsh{"NGSNonRefHitGreater100X"} = "$outDir/flankSeq/NGSNonRefHitGreater100X.fasta";
	$allSeqPathByJunctTypeHsh{"NGSNonRefHitLess1000X"} = "$outDir/flankSeq/NGSNonRefHitLess1000X.fasta";	
	
	$allSeqPathByJunctTypeHsh{"mRNASplicingEffLimA"} = "$outDir/flankSeq/mRNASplicingEffLimA.fasta";
	$allSeqPathByJunctTypeHsh{"mRNASplicingEffLimB"} = "$outDir/flankSeq/mRNASplicingEffLimB.fasta";
	$allSeqPathByJunctTypeHsh{"mRNASplicingEffLimC"} = "$outDir/flankSeq/mRNASplicingEffLimC.fasta";

	$allSeqPathByJunctTypeHsh{"mRNASenseRare"} = "$outDir/flankSeq/mRNASenseRare.fasta";
	$allSeqPathByJunctTypeHsh{"mRNAAntisenseRare"} = "$outDir/flankSeq/mRNAAntisenseRare.fasta";
	$allSeqPathByJunctTypeHsh{"mRNASenseOften"} = "$outDir/flankSeq/mRNASenseOften.fasta";
	$allSeqPathByJunctTypeHsh{"mRNAAntisenseOften"} = "$outDir/flankSeq/mRNAAntisenseOften.fasta";
	$allSeqPathByJunctTypeHsh{"mRNAAntisenseAll"} = "$outDir/flankSeq/mRNAAntisenseAll.fasta";

	$allSeqPathByJunctTypeHsh{"prmntDiff5Seq"} = "$outDir/flankSeq/prmntDiff5Seq.fasta";
	
	open NGSMAJORCP, ">".$allSeqPathByJunctTypeHsh{"NGSMajorClstPrmnt"};
	open NGSMINORNCP, ">".$allSeqPathByJunctTypeHsh{"NGSMinorNonClstPrmnt"};
	open NGSMAJOR, ">".$allSeqPathByJunctTypeHsh{"NGSMajor"};
	open NGSMINOR, ">".$allSeqPathByJunctTypeHsh{"NGSMinor"};
	open NGSCP, ">".$allSeqPathByJunctTypeHsh{"NGSClstPrmnt"};
	open NGSNCP, ">".$allSeqPathByJunctTypeHsh{"NGSNonClstPrmnt"};
	open NGSINFMAJCP, ">".$allSeqPathByJunctTypeHsh{"NGSInfMajorClstPrmnt"};
	open NGSINFMAJNCP, ">".$allSeqPathByJunctTypeHsh{"NGSInfMajorNonClstPrmnt"};

	open NGSINFMAJNCPLIMA, ">".$allSeqPathByJunctTypeHsh{"NGSInfMajorNonClstPrmntLimA"};
	open NGSINFMAJNCPLIMB, ">".$allSeqPathByJunctTypeHsh{"NGSInfMajorNonClstPrmntLimB"};
	open NGSINFMAJNCPLIMC, ">".$allSeqPathByJunctTypeHsh{"NGSInfMajorNonClstPrmntLimC"};
	open NGSINFMAJNCPLIMD, ">".$allSeqPathByJunctTypeHsh{"NGSInfMajorNonClstPrmntLimD"};
	open NGSINFMAJNCPLIME, ">".$allSeqPathByJunctTypeHsh{"NGSInfMajorNonClstPrmntLimE"};

	open NGSSPLCRATIOLIMA, ">".$allSeqPathByJunctTypeHsh{"NGSSplicingRatioLimA"};
	open NGSSPLCRATIOLIMB, ">".$allSeqPathByJunctTypeHsh{"NGSSplicingRatioLimB"};
	open NGSSPLCRATIOLIMC, ">".$allSeqPathByJunctTypeHsh{"NGSSplicingRatioLimC"};
	open NGSSPLCRATIOLIMD, ">".$allSeqPathByJunctTypeHsh{"NGSSplicingRatioLimD"};
	open NGSSPLCRATIOLIME, ">".$allSeqPathByJunctTypeHsh{"NGSSplicingRatioLimE"};

	open NGSSCORELIMA, ">".$allSeqPathByJunctTypeHsh{"NGSScoreLimA"};
	open NGSSCORELIMB, ">".$allSeqPathByJunctTypeHsh{"NGSScoreLimB"};
	open NGSSCORELIMC, ">".$allSeqPathByJunctTypeHsh{"NGSScoreLimC"};
	open NGSSCORELIMD, ">".$allSeqPathByJunctTypeHsh{"NGSScoreLimD"};
	open NGSSCORELIME, ">".$allSeqPathByJunctTypeHsh{"NGSScoreLimE"};
	
	open SPLICEFFLIMA, ">".$allSeqPathByJunctTypeHsh{"mRNASplicingEffLimA"};
	open SPLICEFFLIMB, ">".$allSeqPathByJunctTypeHsh{"mRNASplicingEffLimB"};
	open SPLICEFFLIMC, ">".$allSeqPathByJunctTypeHsh{"mRNASplicingEffLimC"};

	open MRNASR, ">".$allSeqPathByJunctTypeHsh{"mRNASenseRare"};
	open MRNAAR, ">".$allSeqPathByJunctTypeHsh{"mRNAAntisenseRare"};
	open MRNASF, ">".$allSeqPathByJunctTypeHsh{"mRNASenseOften"};
	open MRNAAF, ">".$allSeqPathByJunctTypeHsh{"mRNAAntisenseOften"};
	open MRNAAA, ">".$allSeqPathByJunctTypeHsh{"mRNAAntisenseAll"};

	open NGSREFHIT, ">".$allSeqPathByJunctTypeHsh{"NGSRefHit"};
	open NGSGRTR1000, ">".$allSeqPathByJunctTypeHsh{"NGSNonRefHitGreater1000X"};
	open NGSGRTR100, ">".$allSeqPathByJunctTypeHsh{"NGSNonRefHitGreater100X"};
	open NGSLESS1000, ">".$allSeqPathByJunctTypeHsh{"NGSNonRefHitLess1000X"};
	open NGSLESS1000, ">".$allSeqPathByJunctTypeHsh{"NGSNonRefHitLess1000X"};

	open PRMNTDIFF5SEQ, ">".$allSeqPathByJunctTypeHsh{"prmntDiff5Seq"};

	open PRMNTREADNUM5, ">$outDir/altSpliceSiteShift/prmntReadNumRatio5.txt";
	open PRMNTREADNUM3, ">$outDir/altSpliceSiteShift/prmntReadNumRatio3.txt";
	
	foreach my $junctStr (keys %NGSJunctInfoHsh) {
		
		my $readNum = ${$NGSJunctInfoHsh{$junctStr}}{"readNum"};
		my $score = ${$NGSJunctInfoHsh{$junctStr}}{"score"};
		my $strnd = ${$NGSJunctInfoHsh{$junctStr}}{"strnd"};
		my $leftFlankSeq = ${$NGSJunctInfoHsh{$junctStr}}{"leftFlankSeq"};
		my $rightFlankSeq = ${$NGSJunctInfoHsh{$junctStr}}{"rightFlankSeq"};
		my $upSite = ${$NGSJunctInfoHsh{$junctStr}}{"upSite"};
		my $downSite = ${$NGSJunctInfoHsh{$junctStr}}{"downSite"};
		my $unq = ${$NGSJunctInfoHsh{$junctStr}}{"unq"};
		my $majorJunct = ${$NGSJunctInfoHsh{$junctStr}}{"majorJunct"};
		my $prmntInCluster = ${$NGSJunctInfoHsh{$junctStr}}{"prmntInCluster"};
		my $onPrmntIsofm = ${$NGSJunctInfoHsh{$junctStr}}{"onPrmntIsofm"};
		my $superJunct = ${$NGSJunctInfoHsh{$junctStr}}{"superJunct"};
		my $hitRef = ${$NGSJunctInfoHsh{$junctStr}}{"hitRef"};
		my $cluster = ${$NGSJunctInfoHsh{$junctStr}}{"cluster"};
		my $completeSeq = ${$NGSJunctInfoHsh{$junctStr}}{"completeSeq"};
		my $clstJunctNum = ${$NGSJunctInfoHsh{$junctStr}}{"clstJunctNum"};
		my $prmntReadNumRatio = ${$NGSJunctInfoHsh{$junctStr}}{"prmntReadNumRatio"};
		my $majorCluster = ${$NGSJunctInfoHsh{$junctStr}}{"majorCluster"};
		my $splicingRatio = ${$NGSJunctInfoHsh{$junctStr}}{"splicingRatio"};
		my $BEDCoverRng = ${$NGSJunctInfoHsh{$junctStr}}{"BEDCoverRng"};
		my $ovrlpCtrgy = ${$NGSJunctInfoHsh{$junctStr}}{"ovrlpCtrgy"};
		my $ovrlpSense = ${$NGSJunctInfoHsh{$junctStr}}{"ovrlpSense"};
		my $ovrlpAntisense = ${$NGSJunctInfoHsh{$junctStr}}{"ovrlpAntisense"};
		my $geneCovNt = ${$NGSJunctInfoHsh{$junctStr}}{"geneCovNt"};
		my $senseSplicingEfficiency = ${$NGSJunctInfoHsh{$junctStr}}{"senseSplicingEfficiency"};
		my $antisenseSplicingEfficiency = ${$NGSJunctInfoHsh{$junctStr}}{"antisenseSplicingEfficiency"};
		my $superCtgry = ${$NGSJunctInfoHsh{$junctStr}}{"superCtgry"};
		my $splicingEffType = ${$NGSJunctInfoHsh{$junctStr}}{"splicingEffType"};
		my $prmntDiff5 = ${$NGSJunctInfoHsh{$junctStr}}{"prmntDiff5"};
		my $prmntDiff3 = ${$NGSJunctInfoHsh{$junctStr}}{"prmntDiff3"};
		my $prmntDiffBoth = ${$NGSJunctInfoHsh{$junctStr}}{"prmntDiffBoth"};
		my $consensusValue = ${$NGSJunctInfoHsh{$junctStr}}{"consensusValue"};
		my $intronSize = ${$NGSJunctInfoHsh{$junctStr}}{"intronSize"};
		my $intronFrame = ${$NGSJunctInfoHsh{$junctStr}}{"intronFrame"};

		${${$intronFrameHsh{'NGS'}}{$intronFrame}}{$junctStr} = $intronSize;

		my $refOrAlt = ${$NGSJunctInfoHsh{$junctStr}}{"refOrAlt"};
		my $isofmID = ${$NGSJunctInfoHsh{$junctStr}}{"isofmID"};
		my $altIsofmFullORF = ${$NGSJunctInfoHsh{$junctStr}}{"altIsofmFullORF"};
		my $altSplcType = ${$NGSJunctInfoHsh{$junctStr}}{"altSplcType"};
		my $onRefGene = ${$NGSJunctInfoHsh{$junctStr}}{"onRefGene"};
		my $refJunctAvgReadNum = ${$NGSJunctInfoHsh{$junctStr}}{"refJunctAvgReadNum"};
		my $NGSRefJunctReadNumRatio = ${$NGSJunctInfoHsh{$junctStr}}{"NGSRefJunctReadNumRatio"};
		my $readNumCovPerNtRatio = ${$NGSJunctInfoHsh{$junctStr}}{"readNumCovPerNtRatio"};
		my $altSiteType = ${$NGSJunctInfoHsh{$junctStr}}{"altSiteType"};
		my $exonSkipType = ${$NGSJunctInfoHsh{$junctStr}}{"exonSkipType"};
		my $refShift5 = ${$NGSJunctInfoHsh{$junctStr}}{"refShift5"};
		my $refShift3 = ${$NGSJunctInfoHsh{$junctStr}}{"refShift3"};
		my $skippedExonNum = ${$NGSJunctInfoHsh{$junctStr}}{"skippedExonNum"};
		my $ovrlpRefJunctRdRatio = ${$NGSJunctInfoHsh{$junctStr}}{"ovrlpRefJunctRdRatio"};
		my $avgOvrlpRefJReadNum = ${$NGSJunctInfoHsh{$junctStr}}{"avgOvrlpRefJReadNum"};

		#---Alternative SS
		if (($prmntDiffBoth ne "null") and ($superCtgry eq "mRNA") and ($ovrlpSense ne "null") and ($unq == 1)) {
			${$splicingSiteDiffAllHsh{"prmntDiff5"}}{$prmntDiff5}++;
			${$splicingSiteDiffAllHsh{"prmntDiff3"}}{$prmntDiff3}++;
			${$splicingSiteDiffAllHsh{"prmntDiffBoth"}}{$prmntDiffBoth}++;
			
			print PRMNTREADNUM5 $junctStr."\t".$prmntDiff5."\t".$prmntReadNumRatio."\n" if ($prmntDiff5 != 0);
			print PRMNTREADNUM3 $junctStr."\t".$prmntDiff3."\t".$prmntReadNumRatio."\n" if ($prmntDiff3 != 0);
			
			#---take the prmntDiff5 -4 seq
			my $clusPJ = ${$jClusterInfoHsh{$cluster}}{"prominentJunct"};
			my $PJLeftSeq = ${$NGSJunctInfoHsh{$clusPJ}}{"leftFlankSeq"};
			my $PJRightsSeq = ${$NGSJunctInfoHsh{$clusPJ}}{"rightFlankSeq"};
			my $oriLen = length $PJLeftSeq;
			my $trimLen = $oriLen-$trimExonEnd-$trimIntronEnd;
			die "the trimmed the length of the flank sequence is smaller than zero" if ($trimLen <= 0) ;
			my $trimLeftFlankSeq = substr $PJLeftSeq, $trimExonEnd, $trimLen;
			my $trimRightFlankSeq = substr $PJRightsSeq, $trimIntronEnd, $trimLen;
			my $concateSeq = $trimLeftFlankSeq.$trimRightFlankSeq;

			print PRMNTDIFF5SEQ ">seq\n$concateSeq\n";$junctTypeCountHsh{"prmntDiff5Seq"}++;
		}
		
		#---transfer the info from the prmntJunct to its cluster
		if ($prmntInCluster eq "y") {
		
			${$jClusterInfoHsh{$cluster}}{"ovrlpCtrgy"} = $ovrlpCtrgy;
			${$jClusterInfoHsh{$cluster}}{"ovrlpSense"} = $ovrlpSense; 
			${$jClusterInfoHsh{$cluster}}{"ovrlpAntisense"} = $ovrlpAntisense;  
			${$jClusterInfoHsh{$cluster}}{"superCtgry"} = $superCtgry;
			${$jClusterInfoHsh{$cluster}}{"splicingRatio"} = $splicingRatio;
			${$jClusterInfoHsh{$cluster}}{"splicingEffType"} = $splicingEffType;
		}
		
		print ALLNGSJLOG join "", ((join "\t", ($junctStr, $readNum, $score, $strnd, $unq, $upSite, $downSite, $cluster, $majorCluster, $clstJunctNum, $prmntReadNumRatio, $majorJunct, $prmntInCluster, $onPrmntIsofm, $superJunct, $hitRef, $completeSeq, $leftFlankSeq, $rightFlankSeq, $splicingRatio, $BEDCoverRng, $ovrlpCtrgy, $ovrlpSense, $ovrlpAntisense, $superCtgry, $senseSplicingEfficiency, $antisenseSplicingEfficiency, $prmntDiff5, $prmntDiff3, $prmntDiffBoth, $splicingEffType, $geneCovNt, $refOrAlt, $isofmID, $altIsofmFullORF, $altSplcType, $onRefGene, $refJunctAvgReadNum, $NGSRefJunctReadNumRatio, $readNumCovPerNtRatio, $altSiteType, $exonSkipType, $refShift5, $refShift3, $skippedExonNum, $ovrlpRefJunctRdRatio, $avgOvrlpRefJReadNum, $consensusValue, $intronSize, $intronFrame)), "\n");
		

		#---print the junction seq according to different categoires
		if ($completeSeq eq "complete") {
			my $oriLen = length $leftFlankSeq;
			my $trimLen = $oriLen-$trimExonEnd-$trimIntronEnd;
			die "the trimmed the length of the flank sequence is smaller than zero" if ($trimLen <= 0) ;
			my $trimLeftFlankSeq = substr $leftFlankSeq, $trimExonEnd, $trimLen;
			my $trimRightFlankSeq = substr $rightFlankSeq, $trimIntronEnd, $trimLen;
			my $concateSeq = $trimLeftFlankSeq.$trimRightFlankSeq;
			if (($prmntInCluster eq "y") and ($majorJunct eq "y")) {print NGSMAJORCP ">seq\n$concateSeq\n";$junctTypeCountHsh{"NGSMajorClstPrmnt"}++;push @{$junctTypeScoreHsh{"NGSMajorClstPrmnt"}}, $score;}
			if (($prmntInCluster eq "n") and ($majorJunct eq "n")) {print NGSMINORNCP ">seq\n$concateSeq\n";$junctTypeCountHsh{"NGSMinorNonClstPrmnt"}++;push @{$junctTypeScoreHsh{"NGSMinorNonClstPrmnt"}}, $score;}
			if ($majorJunct eq "y") {print NGSMAJOR ">seq\n$concateSeq\n";$junctTypeCountHsh{"NGSMajor"}++;push @{$junctTypeScoreHsh{"NGSMajor"}}, $score;}
			if ($majorJunct eq "n") {print NGSMINOR ">seq\n$concateSeq\n";$junctTypeCountHsh{"NGSMinor"}++;push @{$junctTypeScoreHsh{"NGSMinor"}}, $score;}
			if (($prmntInCluster eq "y") and ($cluster =~ m/^ic_/)) {print NGSCP ">seq\n$concateSeq\n";$junctTypeCountHsh{"NGSClstPrmnt"}++;push @{$junctTypeScoreHsh{"NGSClstPrmnt"}}, $score;}
			if (($prmntInCluster eq "n") and ($cluster =~ m/^ic_/)) {print NGSNCP ">seq\n$concateSeq\n";$junctTypeCountHsh{"NGSNonClstPrmnt"}++;push @{$junctTypeScoreHsh{"NGSNonClstPrmnt"}}, $score;}
			if (($prmntInCluster eq "y") and ($clstJunctNum > 1) and ($superJunct eq "n") and ($majorCluster eq "y")) {print NGSINFMAJCP ">seq\n$concateSeq\n";$junctTypeCountHsh{"NGSInfMajorClstPrmnt"}++;push @{$junctTypeScoreHsh{"NGSInfMajorClstPrmnt"}}, $score;}
			if (($prmntInCluster eq "n") and ($clstJunctNum > 1) and ($superJunct eq "n") and ($majorCluster eq "y")) {print NGSINFMAJNCP ">seq\n$concateSeq\n";$junctTypeCountHsh{"NGSInfMajorNonClstPrmnt"}++;push @{$junctTypeScoreHsh{"NGSInfMajorNonClstPrmnt"}}, $score;}
			if (($prmntInCluster eq "n") and ($clstJunctNum > 1) and ($superJunct eq "n") and ($majorCluster eq "y") and ($prmntReadNumRatio >= 0.1)) {print NGSINFMAJNCPLIMA ">seq\n$concateSeq\n";$junctTypeCountHsh{"NGSInfMajorNonClstPrmntLimA"}++;push @{$junctTypeScoreHsh{"NGSInfMajorNonClstPrmntLimA"}}, $score;}
			if (($prmntInCluster eq "n") and ($clstJunctNum > 1) and ($superJunct eq "n") and ($majorCluster eq "y") and ($prmntReadNumRatio <= 0.1) and ($prmntReadNumRatio >= 0.01)) {print NGSINFMAJNCPLIMB ">seq\n$concateSeq\n";$junctTypeCountHsh{"NGSInfMajorNonClstPrmntLimB"}++;push @{$junctTypeScoreHsh{"NGSInfMajorNonClstPrmntLimB"}}, $score;}
			if (($prmntInCluster eq "n") and ($clstJunctNum > 1) and ($superJunct eq "n") and ($majorCluster eq "y") and ($prmntReadNumRatio <= 0.01) and ($prmntReadNumRatio >= 0.001)) {print NGSINFMAJNCPLIMC ">seq\n$concateSeq\n";$junctTypeCountHsh{"NGSInfMajorNonClstPrmntLimC"}++;push @{$junctTypeScoreHsh{"NGSInfMajorNonClstPrmntLimC"}}, $score;}
			if (($prmntInCluster eq "n") and ($clstJunctNum > 1) and ($superJunct eq "n") and ($majorCluster eq "y") and ($prmntReadNumRatio <= 0.001) and ($prmntReadNumRatio >= 0.0001)) {print NGSINFMAJNCPLIMD ">seq\n$concateSeq\n";$junctTypeCountHsh{"NGSInfMajorNonClstPrmntLimD"}++;push @{$junctTypeScoreHsh{"NGSInfMajorNonClstPrmntLimD"}}, $score;}
			if (($prmntInCluster eq "n") and ($clstJunctNum > 1) and ($superJunct eq "n") and ($majorCluster eq "y") and ($prmntReadNumRatio <= 0.0001)) {print NGSINFMAJNCPLIME ">seq\n$concateSeq\n";$junctTypeCountHsh{"NGSInfMajorNonClstPrmntLimE"}++;push @{$junctTypeScoreHsh{"NGSInfMajorNonClstPrmntLimE"}}, $score;}

			if ($splicingRatio > 1) {print NGSSPLCRATIOLIMA ">seq\n$concateSeq\n";$junctTypeCountHsh{"NGSSplicingRatioLimA"}++;push @{$junctTypeScoreHsh{"NGSSplicingRatioLimA"}}, $score;}
			if (($splicingRatio <= 1) and ($splicingRatio > 0.1)) {print NGSSPLCRATIOLIMB ">seq\n$concateSeq\n";$junctTypeCountHsh{"NGSSplicingRatioLimB"}++;push @{$junctTypeScoreHsh{"NGSSplicingRatioLimB"}}, $score;}
			if (($splicingRatio <= 0.1) and ($splicingRatio > 0.01)) {print NGSSPLCRATIOLIMC ">seq\n$concateSeq\n";$junctTypeCountHsh{"NGSSplicingRatioLimC"}++;push @{$junctTypeScoreHsh{"NGSSplicingRatioLimC"}}, $score;}
			if (($splicingRatio <= 0.01) and ($splicingRatio > 0.001)) {print NGSSPLCRATIOLIMD ">seq\n$concateSeq\n";$junctTypeCountHsh{"NGSSplicingRatioLimD"}++;push @{$junctTypeScoreHsh{"NGSSplicingRatioLimD"}}, $score;}
			if ($splicingRatio <= 0.001) {print NGSSPLCRATIOLIME ">seq\n$concateSeq\n";$junctTypeCountHsh{"NGSSplicingRatioLimE"}++;push @{$junctTypeScoreHsh{"NGSSplicingRatioLimE"}}, $score;}
			
			if ($score > 1400) {print NGSSCORELIMA ">seq\n$concateSeq\n";$junctTypeCountHsh{"NGSScoreLimA"}++;}
			if (($score <= 1400) and ($score > 1200)) {print NGSSCORELIMB ">seq\n$concateSeq\n";$junctTypeCountHsh{"NGSScoreLimB"}++;}
			if (($score <= 1200) and ($score > 1000)) {print NGSSCORELIMC ">seq\n$concateSeq\n";$junctTypeCountHsh{"NGSScoreLimC"}++;}
			if (($score <= 1000) and ($score > 800)) {print NGSSCORELIMD ">seq\n$concateSeq\n";$junctTypeCountHsh{"NGSScoreLimD"}++;}
			if ($score <= 800) {print NGSSCORELIME ">seq\n$concateSeq\n";$junctTypeCountHsh{"NGSScoreLimE"}++;}

			#my $targetmRNACtgry = "mRNA"; #---for PF
			my $targetmRNACtgry = "bfmRNA"; #---for EHI

			if (($senseSplicingEfficiency ne "null") and ($ovrlpCtrgy eq $targetmRNACtgry)) {
				if ($senseSplicingEfficiency > 0.05) {print SPLICEFFLIMA ">seq\n$concateSeq\n";$junctTypeCountHsh{"mRNASplicingEffLimA"}++;}
				elsif (($senseSplicingEfficiency <= 0.05) and ($senseSplicingEfficiency > 0.005)) {print SPLICEFFLIMB ">seq\n$concateSeq\n";$junctTypeCountHsh{"mRNASplicingEffLimB"}++;}
				elsif ($senseSplicingEfficiency <= 0.005) {print SPLICEFFLIMC ">seq\n$concateSeq\n";$junctTypeCountHsh{"mRNASplicingEffLimC"}++;}
				else {die "Impossible senseSplicingEfficiency value.\n";}
			}

			if ($ovrlpCtrgy eq $targetmRNACtgry) {
				if ($splicingEffType eq "senseRare") {print MRNASR ">seq\n$concateSeq\n";$junctTypeCountHsh{"mRNASenseRare"}++;}
				if ($splicingEffType eq "senseOften") {print MRNASF ">seq\n$concateSeq\n";$junctTypeCountHsh{"mRNASenseOften"}++;}
				if ($splicingEffType eq "antisenseRare") {print MRNAAR ">seq\n$concateSeq\n";$junctTypeCountHsh{"mRNAAntisenseRare"}++;}
				if ($splicingEffType eq "antisenseOften") {print MRNAAF ">seq\n$concateSeq\n";$junctTypeCountHsh{"mRNAAntisenseOften"}++;}
				if (($splicingEffType eq "antisenseRare") or ($splicingEffType eq "antisenseOften")) {print MRNAAA ">seq\n$concateSeq\n";$junctTypeCountHsh{"mRNAAntisenseAll"}++;}
			}
			
			$junctTypeCountHsh{"all"}++; push @{$junctTypeScoreHsh{"all"}}, $score;
			
			if ($hitRef ne "null") {
				print NGSREFHIT ">seq\n$concateSeq\n"; $junctTypeCountHsh{"NGSRefHit"}++;push @{$junctTypeScoreHsh{"NGSRefHit"}}, $score;
			} elsif (($hitRef eq "null") and ($splicingRatio >= 0.01)) {
				print NGSGRTR100 ">seq\n$concateSeq\n"; $junctTypeCountHsh{"NGSNonRefHitGreater100X"}++;push @{$junctTypeScoreHsh{"NGSNonRefHitGreater100X"}}, $score;
			} elsif (($hitRef eq "null") and ($splicingRatio >= 0.001)) {
				print NGSGRTR1000 ">seq\n$concateSeq\n"; $junctTypeCountHsh{"NGSNonRefHitGreater1000X"}++;push @{$junctTypeScoreHsh{"NGSNonRefHitGreater1000X"}}, $score;
			} elsif (($hitRef eq "null") and ($splicingRatio < 0.001)) {
				print NGSLESS1000 ">seq\n$concateSeq\n"; $junctTypeCountHsh{"NGSNonRefHitLess1000X"}++;push @{$junctTypeScoreHsh{"NGSNonRefHitLess1000X"}}, $score;
			} else {
				die;
			}

		}
	}
	
	#---plot
	my $plotCutoff = 50;
	foreach my $bound53 (keys %splicingSiteDiffAllHsh) {
		
		my (%tmpForPlotFullHsh, %tmpForPlotCutoffHsh, @tmpAry);
		
		foreach my $diff (keys  %{$splicingSiteDiffAllHsh{$bound53}}) {
			next if ($diff == 0);
			push @tmpAry, $diff;
			
			$tmpForPlotFullHsh{$diff} = ${$splicingSiteDiffAllHsh{$bound53}}{$diff};

			my $tmpDiff;
			if (($diff < $plotCutoff) and ($diff > -1*$plotCutoff)) {
				$tmpDiff = $diff;
				$tmpForPlotCutoffHsh{$tmpDiff} = ${$splicingSiteDiffAllHsh{$bound53}}{$diff};
			} elsif ($diff >= $plotCutoff) {
				$tmpDiff = $plotCutoff;
				if (not exists $tmpForPlotCutoffHsh{$tmpDiff}) {
					$tmpForPlotCutoffHsh{$tmpDiff} = 0;
				} else {
					$tmpForPlotCutoffHsh{$tmpDiff} += ${$splicingSiteDiffAllHsh{$bound53}}{$diff};
				}
			} elsif ($diff <= -1*$plotCutoff) {
				$tmpDiff = -1*$plotCutoff;
				if (not exists $tmpForPlotCutoffHsh{$tmpDiff}) {
					$tmpForPlotCutoffHsh{$tmpDiff} = 0;
				} else {
					$tmpForPlotCutoffHsh{$tmpDiff} += ${$splicingSiteDiffAllHsh{$bound53}}{$diff};
				}
			} else {
				die;
			}
			
		}
		
		@tmpAry = sort {$a <=> $b} @tmpAry;
		my $min = $tmpAry[0];
		my $max = $tmpAry[-1];
		
		foreach my $value ($min..$max) {
			$tmpForPlotFullHsh{$value} = 0 if (not exists $tmpForPlotFullHsh{$value});
			if (($value <= $plotCutoff) and ($value >= -1*$plotCutoff)) {
				$tmpForPlotCutoffHsh{$value} = 0 if (not exists $tmpForPlotCutoffHsh{$value});
			}
		}
		
		GNUPlotBarChartNumberItem(\%tmpForPlotFullHsh, "intron $bound53 end altSplicingSiteDistance", "$outDir/altSpliceSiteShift/intron$bound53.altSplicingSiteDistance.full.dat", "$outDir/altSpliceSiteShift/intron$bound53.altSplicingSiteDistance.full.pdf", "Distance (nt)", "Frequency", "n", "n");
		GNUPlotBarChartNumberItem(\%tmpForPlotCutoffHsh, "intron $bound53 end altSplicingSiteDistance", "$outDir/altSpliceSiteShift/intron$bound53.altSplicingSiteDistance.cutoff.dat", "$outDir/altSpliceSiteShift/intron$bound53.altSplicingSiteDistance.cutoff.pdf", "Distance (nt)", "Frequency", "n", "n");
	}
	

	#---print all info of the intronClusters and Junctions
	my %mRNASenseClusterCountHsh;
	open (ALLCLSTINFO, ">$outDir/junctInfo/allJunctionClusterInfo.txt");
	open (readNumALTJTOTAL, ">$outDir/correlation/readNumVsAltJunctInClusterTotal.txt"); 
	open (readNumALTJMAJOR, ">$outDir/correlation/readNumVsAltJunctInClusterMajor.txt");
	open (readNumCONVSALT, ">$outDir/correlation/readNumAltJunctVsConstitutiveJunct.txt");
	
	print ALLCLSTINFO join "", ((join "\t", ("clusterName", "superCluster", "hitRef", "unq", "majorJunctNum", "totalJunctNum", "majorCluster", "prominentJunct", "mostSupportReadNum", "totalRdNumInCluster", "allJunctStr", "ovrlpAS", "ovrlpSS", "ovrlpCtrgy", "ovrlpSense", "ovrlpAntisense", "superCtgry", "splicingRatio")), "\n");
	foreach my $clusterName (sort {$a cmp $b} keys %jClusterInfoHsh) {

		my $ovrlpCtrgy = ${$jClusterInfoHsh{$clusterName}}{"ovrlpCtrgy"};
		my $ovrlpSense = ${$jClusterInfoHsh{$clusterName}}{"ovrlpSense"};
		my $ovrlpAntisense = ${$jClusterInfoHsh{$clusterName}}{"ovrlpAntisense"};
		my $superCtgry = ${$jClusterInfoHsh{$clusterName}}{"superCtgry"};
		my $splicingRatio = ${$jClusterInfoHsh{$clusterName}}{"splicingRatio"};
		my $splicingEffType = ${$jClusterInfoHsh{$clusterName}}{"splicingEffType"};
		my $hitRef = ${$jClusterInfoHsh{$clusterName}}{"hitRef"};
		my $unq = ${$jClusterInfoHsh{$clusterName}}{"unq"};
		my $majorJunctNum = ${$jClusterInfoHsh{$clusterName}}{"majorJunctNum"};
		my $totalJunctNum = ${$jClusterInfoHsh{$clusterName}}{"totalJunctNum"};
		my $majorCluster = ${$jClusterInfoHsh{$clusterName}}{"majorCluster"};
		my $prominentJunct = ${$jClusterInfoHsh{$clusterName}}{"prominentJunct"};
		my $mostSupportReadNum = ${$jClusterInfoHsh{$clusterName}}{"mostSupportReadNum"};
		my $allJunctStr = join ";", @{${$jClusterInfoHsh{$clusterName}}{"allJunctStr"}};
		my $ovrlpAS = join ";", @{${$jClusterInfoHsh{$clusterName}}{"ovrlpAS"}};
		my $ovrlpSS = join ";", @{${$jClusterInfoHsh{$clusterName}}{"ovrlpSS"}};
		my $superCluster = ${$jClusterInfoHsh{$clusterName}}{"superCluster"};
		my $totalRdNumInCluster = ${$jClusterInfoHsh{$clusterName}}{"totalRdNumInCluster"};
		
		if (($superCtgry eq "mRNA") and ($ovrlpSense ne "null")) {#---overlapping with mRNA and sense
			
			$mRNASenseClusterCountHsh{"total"}++;
			${$mRNASenseClusterCountHsh{"majorCluster"}}{$majorCluster}++;
			${$mRNASenseClusterCountHsh{"totalJunctNum"}}{$totalJunctNum}++;
			push @{${$mRNASenseClusterCountHsh{"readNumAryByJunctNum"}}{$totalJunctNum}}, $mostSupportReadNum;
			print readNumALTJTOTAL $totalJunctNum."\t".$mostSupportReadNum."\n" if ($mostSupportReadNum >= 20);
			if ($majorCluster eq "y") {
				print readNumALTJMAJOR $totalJunctNum."\t".$mostSupportReadNum."\n" if ($mostSupportReadNum >= 20);
				${$mRNASenseClusterCountHsh{"totalJunctNumMajor"}}{$totalJunctNum}++;
				push @{${$mRNASenseClusterCountHsh{"readNumAryByJunctNumMajor"}}{$totalJunctNum}}, $mostSupportReadNum;
				
				#----get the constitutive and alternative junction num of read
				my $altReadNum = my $constitutiveReadNum = 0;
				my $altJunctNum = @{${$jClusterInfoHsh{$clusterName}}{"allJunctStr"}} - 1;
				foreach my $junctStr (@{${$jClusterInfoHsh{$clusterName}}{"allJunctStr"}}) {
					my $readNum = ${$NGSJunctInfoHsh{$junctStr}}{"readNum"};
					if ($junctStr eq $prominentJunct) {
						$constitutiveReadNum = $readNum;
					} else {
						$altReadNum += $readNum;
					}
				}
				print readNumCONVSALT $constitutiveReadNum."\t".$altReadNum."\t".$altJunctNum."\n";
			}
		}
		
		print ALLCLSTINFO join "", ((join "\t", ($clusterName, $superCluster, $hitRef, $unq, $majorJunctNum, $totalJunctNum, $majorCluster, $prominentJunct, $mostSupportReadNum, $totalRdNumInCluster, $allJunctStr, $ovrlpAS, $ovrlpSS, $ovrlpCtrgy, $ovrlpSense, $ovrlpAntisense, $superCtgry, $splicingRatio)), "\n");
	}
	close ALLCLSTINFO;
	
	print "\nCounting of mRNA sense junction clusters\n";
	print "Total=".$mRNASenseClusterCountHsh{"total"}."\tMajor=".${$mRNASenseClusterCountHsh{"majorCluster"}}{"y"}."\tMinor=".${$mRNASenseClusterCountHsh{"majorCluster"}}{"n"}."\n";
	print "JunctNumInAllClust\tNum\tPct\n";
	foreach my $totalJunctNum (sort {$a <=> $b} keys %{$mRNASenseClusterCountHsh{"totalJunctNum"}}) {
		my $count = ${$mRNASenseClusterCountHsh{"totalJunctNum"}}{$totalJunctNum};
		my $pct = sprintf "%0.2f", 100*$count/$mRNASenseClusterCountHsh{"total"};
		my ($mean, $SD) = calculateStandardDeviationAndMean(\@{${$mRNASenseClusterCountHsh{"readNumAryByJunctNum"}}{$totalJunctNum}});
		print $totalJunctNum."\t".$count."\t".$pct."\t".$mean."\t".$SD."\n";
	}

	print "JunctNumInMajorClust\tNum\tPct\n";
	foreach my $totalJunctNum (sort {$a <=> $b} keys %{$mRNASenseClusterCountHsh{"totalJunctNumMajor"}}) {
		my $count = ${$mRNASenseClusterCountHsh{"totalJunctNumMajor"}}{$totalJunctNum};
		my $pct = sprintf "%0.2f", 100*$count/${$mRNASenseClusterCountHsh{"majorCluster"}}{"y"};
		my ($mean, $SD) = calculateStandardDeviationAndMean(\@{${$mRNASenseClusterCountHsh{"readNumAryByJunctNumMajor"}}{$totalJunctNum}});
		print $totalJunctNum."\t".$count."\t".$pct."\t".$mean."\t".$SD."\n";
	}
	
	close ALLNGSJLOG;
	close NGSMAJORCP;
	close NGSMINORNCP;
	close NGSMAJOR;
	close NGSMINOR;
	close NGSCP;
	close NGSNCP;
	close NGSINFMAJCP;
	close NGSINFMAJNCP;
	close NGSINFMAJNCPLIMA;
	close NGSINFMAJNCPLIMB;
	close NGSINFMAJNCPLIMC;
	close NGSINFMAJNCPLIMD;
	close NGSINFMAJNCPLIME;
	close NGSSPLCRATIOLIMA;
	close NGSSPLCRATIOLIMB;
	close NGSSPLCRATIOLIMC;
	close NGSSPLCRATIOLIMD;
	close NGSSPLCRATIOLIME;
	

	open ALLREFJLOG, ">$outDir/junctInfo/all_REF_junct_log.txt";
	print ALLREFJLOG join "", ((join "\t", ("junctStr", "NGSHit", "comment", "unq", "strnd", "upSite", "downSite", "completeSeq", "leftFlankSeq", "rightFlankSeq", "confirmedNGSJReadNum", "ovlpNGSJInferior", "refGeneCovPerNt", "readNumCovPerNtRatio", "intronSize", "intronFrame")), "\n");
	my %refJuncthitSceneCountHsh;
	foreach my $junctStr (keys %refJunctInfoHsh) {
		
		my $NGSHit = ${$refJunctInfoHsh{$junctStr}}{"NGSHit"};
		my $unq = ${$refJunctInfoHsh{$junctStr}}{"unq"};
		my $strnd = ${$refJunctInfoHsh{$junctStr}}{"strnd"};
		my $leftFlankSeq = ${$refJunctInfoHsh{$junctStr}}{"leftFlankSeq"};
		my $rightFlankSeq = ${$refJunctInfoHsh{$junctStr}}{"rightFlankSeq"};
		my $upSite = ${$refJunctInfoHsh{$junctStr}}{"upSite"};
		my $downSite = ${$refJunctInfoHsh{$junctStr}}{"downSite"};
		my $completeSeq = ${$refJunctInfoHsh{$junctStr}}{"completeSeq"};
		my $intronSize = ${$refJunctInfoHsh{$junctStr}}{"intronSize"};
		my $intronFrame = ${$refJunctInfoHsh{$junctStr}}{"intronFrame"};
		
		${${$intronFrameHsh{'ref'}}{$intronFrame}}{$junctStr} = $intronSize;
		
		my $comment = "null";

		my $confirmedNGSJReadNum = ${$refJunctInfoHsh{$junctStr}}{"confirmedNGSJReadNum"};
		my $ovlpNGSJInferior = ${$refJunctInfoHsh{$junctStr}}{"ovlpNGSJInferior"};
		my $refGeneCovPerNt = ${$refJunctInfoHsh{$junctStr}}{"refGeneCovPerNt"};
		my $readNumCovPerNtRatio = ${$refJunctInfoHsh{$junctStr}}{"readNumCovPerNtRatio"};
		
		$comment = $refJunctCommentHsh{$junctStr} if (exists $refJunctCommentHsh{$junctStr});
		$refJuncthitSceneCountHsh{$NGSHit}++;
		print ALLREFJLOG join "", ((join "\t", ($junctStr, $NGSHit, $comment, $unq, $strnd, $upSite, $downSite, $completeSeq, $leftFlankSeq, $rightFlankSeq, $confirmedNGSJReadNum, $ovlpNGSJInferior, $refGeneCovPerNt, $readNumCovPerNtRatio, $intronSize, $intronFrame)), "\n");
	}
	close ALLREFJLOG;


	print "print the intronFrameCount\n";
	open INTRONFRAMELOG, ">$outDir/junctInfo/intron.3n.frame.log.txt";
	foreach my $refOrNGS (sort {$a cmp $b} keys %intronFrameHsh) {
		foreach my $intronFrame (sort {$a cmp $b} keys %{$intronFrameHsh{$refOrNGS}}) {
			my $junctCount = keys %{${$intronFrameHsh{$refOrNGS}}{$intronFrame}};
			print INTRONFRAMELOG join "", ((join "\t", ($refOrNGS, $intronFrame, $junctCount)), "\n");
		}
	}
	close INTRONFRAMELOG;

	print "Final count of Ref-NGS junction Hit\n";
	my $totalRefJunctNum = keys %refJunctInfoHsh;
	foreach my $scene (sort {$refJuncthitSceneCountHsh{$b} <=> $refJuncthitSceneCountHsh{$a}} keys %refJuncthitSceneCountHsh) {
		my $pct = sprintf "%.06f", 100*$refJuncthitSceneCountHsh{$scene}/$totalRefJunctNum;
		print $scene."\t".$refJuncthitSceneCountHsh{$scene}."\t".$pct."\n";
	}
	
	print "Final count of NGS junction types\n";
	foreach my $junctType (sort {$a cmp $b} keys %junctTypeCountHsh) {
		print $junctType."\t".$junctTypeCountHsh{$junctType}."\n";
	}
	
	#----generate the score histogram data file
	my @tmpAry;
	foreach my $junctType (sort {$a cmp $b} keys %junctTypeScoreHsh) {
		push @tmpAry, @{$junctTypeScoreHsh{$junctType}};
	}	
	
	@tmpAry = sort {$a <=> $b} @tmpAry;
	my $min = int $tmpAry[0];
	my $max = int $tmpAry[-1];
	my $bin_width = int (($max-$min)/100);

	my (%allScoreFreqByBinHsh, %allScorePrptnByBinHsh, @junctTypeAry);
	foreach my $junctType (sort {$a cmp $b} keys %junctTypeScoreHsh) {
		open SCOREFREQPROP, ">$outDir/scoreReadNum/$junctType.Score.Freq.Proportion.txt";
		my $tmpFreqByBinHsh_ref = histogram($bin_width, $max, $min, \@{$junctTypeScoreHsh{$junctType}});
		my %tmpFreqByBinHsh = %{$tmpFreqByBinHsh_ref};
		push @junctTypeAry, $junctType;
		
		foreach my $bin (sort {$a <=> $b} keys %tmpFreqByBinHsh) {
			my $freq = $tmpFreqByBinHsh{$bin};
			my $totalNum = $junctTypeCountHsh{$junctType};
			my $proportion = sprintf "%0.5f", ($freq/$totalNum);
			print SCOREFREQPROP $bin."\t".$freq."\t".$proportion."\n";
			push @{$allScoreFreqByBinHsh{$bin}}, $freq;
			push @{$allScorePrptnByBinHsh{$bin}}, $proportion;
		}
		close SCOREFREQPROP;
	}	
	
	open SCOREFREQ, ">$outDir/scoreReadNum/all.Freq.PrismReady.txt";
	open SCOREPROP, ">$outDir/scoreReadNum/all.Proportion.PrismReady.txt";
	print SCOREFREQ join "", ((join "\t", ("bin", @junctTypeAry)), "\n");
	print SCOREPROP join "", ((join "\t", ("bin", @junctTypeAry)), "\n");

	foreach my $bin (sort {$a <=> $b} keys %allScoreFreqByBinHsh) {
		print SCOREFREQ join "", ((join "\t", ($bin, @{$allScoreFreqByBinHsh{$bin}})), "\n");
		print SCOREPROP join "", ((join "\t", ($bin, @{$allScorePrptnByBinHsh{$bin}})), "\n");
	}
	
	close SCOREPROP;
	close SCOREFREQ;
	
	return (\%allSeqPathByJunctTypeHsh);
	
}
########################################################################## readRefJunctComment
sub readRefJunctComment {

	my %refJunctCommentHsh;
	if ($missModJunctComPath ne "no") {
		open INFILE, "$missModJunctComPath";
		while (my $theLine = <INFILE>) {
			chomp $theLine;
			my @theLineSplt = split /\t/, $theLine;
			my $junctStr = $theLineSplt[0];
			my $comment = $theLineSplt[2];
			$refJunctCommentHsh{$junctStr} = $comment;
		}
		close INFILE;
	}
	
	return \%refJunctCommentHsh;
	
}
########################################################################## calculateEntropyAndPlotLogo
sub calculateEntropyAndPlotLogo {
	
	my %allSeqPathByJunctTypeHsh = %{$_[0]};
	my $length = 50;
	
	my $junctTypeNum = keys %allSeqPathByJunctTypeHsh;
	
	my $progCount=0;
	my $allPdfPath = "";
	my $errorLogPath = "$outDir/flankSeq/error.log.txt";
	open (ERROR, ">$errorLogPath");
	close ERROR;
	my $filterSeqLogPath = "$outDir/flankSeq/filtered.JunctInfo.txt";

	foreach my $junctSeqType (sort {$b cmp $a} keys %allSeqPathByJunctTypeHsh) {
		$progCount++;

		my $entropyPath = "$outDir/flankSeq/$junctSeqType.entropy.txt";
		my $seqPath = $allSeqPathByJunctTypeHsh{$junctSeqType};
		my $logoPath = "$outDir/flankSeq/$junctSeqType.weblogo.pdf";
		
		$allPdfPath .= " \'$logoPath\'";
		
		print "Running Weblogo for $junctSeqType.\n";
		
		my $title = $junctSeqType;
		my $numSeq = `sed -n \'\$=\' $seqPath`;
		$numSeq = $numSeq/2;
		$title .= "_n=$numSeq";
		#my $weblogoEntropyCMD = "weblogo -f $seqPath -F txt -n $length -s large --title $title -A rna -c classic >$entropyPath";
		my $weblogoEntropyCMD = "weblogo -f $seqPath -F logodata -n $length -s large --title $title -A rna -c classic >$entropyPath";
		my $grepCMD = "ps -ef | grep weblogo | grep -v grep | grep -v perl";
		runAndCheckSerialTask($grepCMD, "weblogo", $weblogoEntropyCMD, "$outDir/error.log.txt");

		#$weblogoEntropyCMD = "weblogo -f $seqPath -F pdf -n $length -s large --title $title -A rna -c classic >$logoPath";
		$weblogoEntropyCMD = "weblogo -f $seqPath -F pdf -n $length -s large --title $title -A rna -c classic >$logoPath";
		$grepCMD = "ps -ef | grep weblogo | grep -v grep | grep -v perl";
		runAndCheckSerialTask($grepCMD, "weblogo", $weblogoEntropyCMD, $errorLogPath);

	}
	
	print "Merging all pdfs.\n";
	my $mergePdfCMD = "gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=\'$outDir/flankSeq/all.merged.pdf\' $allPdfPath";
	my $grepCMD = "ps -ef | grep $allPdfPath | grep -v grep";
	runAndCheckSerialTask($grepCMD, $mergePdfCMD, $mergePdfCMD, $errorLogPath);
	
}
########################################################################## runAndCheckSerialTask
sub runAndCheckSerialTask {

	my $grepCmd = $_[0];
	my $grepStr = $_[1];
	my $cmd = $_[2];
	my $errorLogPath = $_[3];

	system (qq|$cmd 2>>$errorLogPath &|);
	my $sdout = $grepStr;
	while ($sdout =~ m/$grepStr/) {
		$sdout = `$grepCmd`;
		sleep (0.001);
	}
}
########################################################################## histogram
sub histogram {

   	my $bin_width = $_[0];
   	my $max = $_[1];
   	my $min = $_[2];
   	my @list = @{$_[3]};
	
	my %freqByBinHsh;
	
	for (my $bin = $min; $bin <= $max; $bin += $bin_width) {
		my $freq = 0;
		foreach my $value (@list) {
			$freq++ if (($value <= $bin+$bin_width) and ($value > $bin));
		}
		$freqByBinHsh{$bin} = $freq;
	}	

	return \%freqByBinHsh;
 }
########################################################################## countNGSJunctOverlapFtur
sub countNGSJunctOverlapFtur {
	
	my %NGSJunctInfoHsh = %{$_[0]};
	my %SSOvrlpNGSJFturCountHsh = %{$_[1]};
	my %XSOvrlpNGSJFturCountHsh = %{$_[2]};
	my %fturCountCtgryByGeneHsh = %{$_[3]};
	my %countFturCovInfoHsh = %{$_[4]};
	my %junctReadNumHsh = %{$_[5]};
	my %SSOvrlpJunctIntronHsh = %{$_[6]};
	
	my (%ctgryHitCountHsh, %ctgryCountHsh, %superCtgryGeneCountHsh, %superCtgryJunctCountHsh, %ctgrySenseRegularSpliceGeneCountHsh, %superCtgrySenseRegularSpliceGeneCountHsh, @hitRefmRNASenseSplicingEffHshAry, %geneBasedJunctCountHsh);
	
	#---define the superCtgry
	my %superCtgryDefineHsh;
	
	#---EHI categories
	$superCtgryDefineHsh{"rRNA"}= "rRNA";
	$superCtgryDefineHsh{"reprDNA"}= "reprDNA";
	$superCtgryDefineHsh{"tRNA"}= "tRNA";
	$superCtgryDefineHsh{"uTrnsfg"}= "ncRNA";
	$superCtgryDefineHsh{"dTrnsfg"}= "ncRNA";
	$superCtgryDefineHsh{"orimRNA"}= "mRNA";
	$superCtgryDefineHsh{"bfdnmRNA"}= "mRNA";
	$superCtgryDefineHsh{"bfdomRNA"}= "mRNA";
	$superCtgryDefineHsh{"bfunmRNA"}= "mRNA";
	$superCtgryDefineHsh{"bfuomRNA"}= "mRNA";
	$superCtgryDefineHsh{"bfmRNA"}= "mRNA";
	$superCtgryDefineHsh{"ERE"}= "repElmnt";
	$superCtgryDefineHsh{"EhRC"}= "repElmnt";
	$superCtgryDefineHsh{"LINE"}= "repElmnt";
	$superCtgryDefineHsh{"SINE"}= "repElmnt";
	$superCtgryDefineHsh{"null"}= "null";

=pod
	#---PF categories
	$superCtgryDefineHsh{"rRNA"}= "rRNA";
	$superCtgryDefineHsh{"plstd_rRNA"}= "rRNA";
	$superCtgryDefineHsh{"mito_rRNA"}= "rRNA";
	$superCtgryDefineHsh{"tRNA"}= "tRNA";
	$superCtgryDefineHsh{"plstd_tRNA"}= "tRNA";
	$superCtgryDefineHsh{"ncRNA"}= "ncRNA";
	$superCtgryDefineHsh{"snoRNA"}= "snoRNA";
	$superCtgryDefineHsh{"snRNA"}= "snRNA";
	$superCtgryDefineHsh{"pseudogenic_transcript"}= "pseudogenic_transcript";
	$superCtgryDefineHsh{"mRNA"}= "mRNA";
	$superCtgryDefineHsh{"plstd_mRNA"}= "mRNA";
	$superCtgryDefineHsh{"mito_mRNA"}= "mRNA";
	$superCtgryDefineHsh{"null"}= "null";
=cut

	foreach my $ctgry (keys %superCtgryDefineHsh) {
		my $superCtgry = $superCtgryDefineHsh{$ctgry};
		${$superCtgryGeneCountHsh{$superCtgry}}{"total"} = 0;
		${$superCtgryGeneCountHsh{$superCtgry}}{"sense"} = 0;
		${$superCtgryGeneCountHsh{$superCtgry}}{"antisense"} = 0;
		%{${$superCtgrySenseRegularSpliceGeneCountHsh{$superCtgry}}{"rare"}} = ();
		%{${$superCtgrySenseRegularSpliceGeneCountHsh{$superCtgry}}{"often"}} = ();
	}
	
	foreach my $geneID (keys %fturCountCtgryByGeneHsh) {
		my $ctgry = $fturCountCtgryByGeneHsh{$geneID};
		${${$ctgryHitCountHsh{$ctgry}}{"total"}}{$geneID}++;
		%{${$ctgryHitCountHsh{$ctgry}}{"sense"}} = ();
		%{${$ctgryHitCountHsh{$ctgry}}{"sense"}} = ();
		%{${$ctgrySenseRegularSpliceGeneCountHsh{$ctgry}}{"rare"}} = ();
		%{${$ctgrySenseRegularSpliceGeneCountHsh{$ctgry}}{"often"}} = ();
	}
	
	#---define what is "rare" what is "often" based on the 95% percentile of all junct with hitRef on bfuomRNA and bfunmRNA
	foreach my $junctStr (keys %SSOvrlpNGSJFturCountHsh) {
		foreach my $geneID (keys %{$SSOvrlpNGSJFturCountHsh{$junctStr}}) {
			my $ctgry = $fturCountCtgryByGeneHsh{$geneID};
			my $readNum = $junctReadNumHsh{$junctStr};
			my $geneCovNt = (${$countFturCovInfoHsh{$geneID}}{"plusCov"} + ${$countFturCovInfoHsh{$geneID}}{"minusCov"}) / ${$countFturCovInfoHsh{$geneID}}{"length"};
			my $splicingEfficiency = 999;
			$splicingEfficiency = $readNum/$geneCovNt if ($geneCovNt > 0);
			if ((${$NGSJunctInfoHsh{$junctStr}}{"hitRef"} ne "null") and (($ctgry eq "mRNA") or ($ctgry eq "bfmRNA") or ($ctgry eq "bfumRNA"))) {#----mRNA fpr PF and others, bfmRNA and bfumRNA for EHI
				push @hitRefmRNASenseSplicingEffHshAry, $splicingEfficiency;
			}
		}	
	}
	
	@hitRefmRNASenseSplicingEffHshAry = sort {$b <=> $a} @hitRefmRNASenseSplicingEffHshAry;
	
	my $totalBfumRNANum = @hitRefmRNASenseSplicingEffHshAry;
	my $pctCut = int ($pctCutLim*$totalBfumRNANum);
	my $oftenSplicingEffCutoff = $hitRefmRNASenseSplicingEffHshAry[$pctCut];
	
	print "$oftenSplicingEffCutoff is the $pctCutLim percentile of all junct with hitRef on mRNA, used as oftenSplicingEffCutoff\n";
	
	my $totalJunctNum = keys %NGSJunctInfoHsh;
	
	#---go through all juncts
	foreach my $junctStr (keys %NGSJunctInfoHsh) {
		my $readNum = $junctReadNumHsh{$junctStr};
		${$NGSJunctInfoHsh{$junctStr}}{"senseSplicingEfficiency"} = "null";
		${$NGSJunctInfoHsh{$junctStr}}{"antisenseSplicingEfficiency"} = "null";
		${$NGSJunctInfoHsh{$junctStr}}{"splicingEffType"} = "null";
		${$NGSJunctInfoHsh{$junctStr}}{"geneCovNt"} = "null";
		my $splicingEfficiency = -1;
		my @ctrgyAry = my @senseAry = my @antisenseAry;
		my $splicingEffType = "null";

		if (exists 	$XSOvrlpNGSJFturCountHsh{$junctStr}) {#---overlap with something
			my $geneNum = keys %{$XSOvrlpNGSJFturCountHsh{$junctStr}};
			foreach my $geneID (keys %{$XSOvrlpNGSJFturCountHsh{$junctStr}}) {
				my $ctgry = $fturCountCtgryByGeneHsh{$geneID};
				my $superCtgry = $superCtgryDefineHsh{$ctgry};
				if ($geneNum == 1) {#---record junctions only overlap with one ftur
					my $geneCovNt = (${$countFturCovInfoHsh{$geneID}}{"plusCov"} + ${$countFturCovInfoHsh{$geneID}}{"minusCov"}) / ${$countFturCovInfoHsh{$geneID}}{"length"};
					my $splicingEfficiency = 999;
					$splicingEfficiency = $readNum/$geneCovNt if ($geneCovNt > 0);
					${$NGSJunctInfoHsh{$junctStr}}{"geneCovNt"} = $geneCovNt;
					
					if (exists ${$SSOvrlpNGSJFturCountHsh{$junctStr}}{$geneID}) {
						${$NGSJunctInfoHsh{$junctStr}}{"senseSplicingEfficiency"} = $splicingEfficiency;
						push @{${${$ctgryHitCountHsh{$ctgry}}{"sense"}}{$geneID}}, $splicingEfficiency;

						if ($splicingEfficiency < $oftenSplicingEffCutoff) {
							$splicingEffType = "senseRare";
							
						} else {#----> greater than $oftenSplicingEffCutoff
							$splicingEffType = "senseOften";
						}
				
						#----collect gene based info
						${$geneBasedJunctCountHsh{$geneID}}{"senseTotal"}++;
						${$geneBasedJunctCountHsh{$geneID}}{$splicingEffType}++;

					} else {
						
						if ($splicingEfficiency < $oftenSplicingEffCutoff) {
							$splicingEffType = "antisenseRare";
						} else {#----> greater than $oftenSplicingEffCutoff
							$splicingEffType = "antisenseOften";
						}
						
						#----collect gene based info
						${$geneBasedJunctCountHsh{$geneID}}{"antisenseTotal"}++;
						${$geneBasedJunctCountHsh{$geneID}}{$splicingEffType}++;

						${$NGSJunctInfoHsh{$junctStr}}{"antisenseSplicingEfficiency"} = $splicingEfficiency;
						push @{${${$ctgryHitCountHsh{$ctgry}}{"antisense"}}{$geneID}}, $splicingEfficiency;
					}
					

					${${$ctgrySenseRegularSpliceGeneCountHsh{$ctgry}}{$splicingEffType}}{$geneID}++;
					${${$superCtgrySenseRegularSpliceGeneCountHsh{$superCtgry}}{$splicingEffType}}{$geneID}++;
					${$NGSJunctInfoHsh{$junctStr}}{"splicingEffType"} = $splicingEffType;
				}
				
				push @ctrgyAry, $ctgry;
				
				if (exists ${$SSOvrlpNGSJFturCountHsh{$junctStr}}{$geneID}) {#---sense
					push @senseAry, $geneID;
				} else {#---antisense
					push @antisenseAry, $geneID;
				}
			}
		}
		
		push @ctrgyAry, "null" if (@ctrgyAry < 1);
		push @senseAry, "null" if (@senseAry < 1);
		push @antisenseAry, "null" if (@antisenseAry < 1);
		
		${$NGSJunctInfoHsh{$junctStr}}{"ovrlpCtrgy"} = join ";", @ctrgyAry;
		${$NGSJunctInfoHsh{$junctStr}}{"ovrlpSense"} = join ";", @senseAry;
		${$NGSJunctInfoHsh{$junctStr}}{"ovrlpAntisense"} = join ";", @antisenseAry;
		my $ctgry = ${$NGSJunctInfoHsh{$junctStr}}{"ovrlpCtrgy"};
		${$ctgryCountHsh{$ctgry}}{"total"}++;
		${$ctgryCountHsh{$ctgry}}{"sense"} = 0 if (not exists ${$ctgryCountHsh{$ctgry}}{"sense"});
		${$ctgryCountHsh{$ctgry}}{"antisense"} = 0 if (not exists ${$ctgryCountHsh{$ctgry}}{"antisense"});
		${$ctgryCountHsh{$ctgry}}{"senseOften"} = 0 if (not exists ${$ctgryCountHsh{$ctgry}}{"senseOften"});
		${$ctgryCountHsh{$ctgry}}{"senseRare"} = 0 if (not exists ${$ctgryCountHsh{$ctgry}}{"senseRare"});
		${$ctgryCountHsh{$ctgry}}{"antisenseOften"} = 0 if (not exists ${$ctgryCountHsh{$ctgry}}{"antisenseOften"});
		${$ctgryCountHsh{$ctgry}}{"antisenseRare"} = 0 if (not exists ${$ctgryCountHsh{$ctgry}}{"antisenseRare"});
		${$ctgryCountHsh{$ctgry}}{"null"} = 0 if (not exists ${$ctgryCountHsh{$ctgry}}{"null"});
		
		my $superCtgry = "moreThanOne";
		$superCtgry = $superCtgryDefineHsh{$ctgry} if (exists $superCtgryDefineHsh{$ctgry});
		${$NGSJunctInfoHsh{$junctStr}}{"superCtgry"} = $superCtgry;

		${$superCtgryJunctCountHsh{$superCtgry}}{"total"}++;
		${$superCtgryJunctCountHsh{$superCtgry}}{"sense"} = 0 if (not exists ${$superCtgryJunctCountHsh{$superCtgry}}{"sense"});
		${$superCtgryJunctCountHsh{$superCtgry}}{"antisense"} = 0 if (not exists ${$superCtgryJunctCountHsh{$superCtgry}}{"antisense"});
		${$superCtgryJunctCountHsh{$superCtgry}}{"senseRare"} = 0 if (not exists ${$superCtgryJunctCountHsh{$superCtgry}}{"senseRare"});
		${$superCtgryJunctCountHsh{$superCtgry}}{"senseOften"} = 0 if (not exists ${$superCtgryJunctCountHsh{$superCtgry}}{"senseOften"});
		${$superCtgryJunctCountHsh{$superCtgry}}{"antisenseRare"} = 0 if (not exists ${$superCtgryJunctCountHsh{$superCtgry}}{"antisenseRare"});
		${$superCtgryJunctCountHsh{$superCtgry}}{"antisenseOften"} = 0 if (not exists ${$superCtgryJunctCountHsh{$superCtgry}}{"antisenseOften"});
		${$superCtgryJunctCountHsh{$superCtgry}}{"null"} = 0 if (not exists ${$superCtgryJunctCountHsh{$superCtgry}}{"null"});
		
		if (${$NGSJunctInfoHsh{$junctStr}}{"ovrlpSense"} ne "null") {
			${$ctgryCountHsh{$ctgry}}{"sense"}++;
			${$superCtgryJunctCountHsh{$superCtgry}}{"sense"}++;
		}

		if (${$NGSJunctInfoHsh{$junctStr}}{"ovrlpAntisense"} ne "null") {
			${$ctgryCountHsh{$ctgry}}{"antisense"}++;
			${$superCtgryJunctCountHsh{$superCtgry}}{"antisense"}++;
		}

		${$ctgryCountHsh{$ctgry}}{$splicingEffType}++;
		${$superCtgryJunctCountHsh{$superCtgry}}{$splicingEffType}++;
		
	}
	
	#---check mRNA rare freqent overlapping
	my %rareOftenOvrlpCountHsh;
	foreach my $junctStr (keys %NGSJunctInfoHsh) {
		my @tmpAry = ();
		my $ctgry = ${$NGSJunctInfoHsh{$junctStr}}{"ovrlpCtrgy"};
		my $superCtgry = "moreThanOne";
		$superCtgry = $superCtgryDefineHsh{$ctgry} if (exists $superCtgryDefineHsh{$ctgry});
		${$rareOftenOvrlpCountHsh{$superCtgry}}{"senseRare"} = 0 if (not exists ${$rareOftenOvrlpCountHsh{$superCtgry}}{"senseRare"});
		${$rareOftenOvrlpCountHsh{$superCtgry}}{"senseOften"} = 0 if (not exists ${$rareOftenOvrlpCountHsh{$superCtgry}}{"senseOften"});
		${$rareOftenOvrlpCountHsh{$superCtgry}}{"antisenseRare"} = 0 if (not exists ${$rareOftenOvrlpCountHsh{$superCtgry}}{"antisenseRare"});
		${$rareOftenOvrlpCountHsh{$superCtgry}}{"antisenseOften"} = 0 if (not exists ${$rareOftenOvrlpCountHsh{$superCtgry}}{"antisenseOften"});

		my $splicingEffType = ${$NGSJunctInfoHsh{$junctStr}}{"splicingEffType"}; 
		if (exists $SSOvrlpJunctIntronHsh{$junctStr}) {
			foreach my $hitJunctStr (keys %{$SSOvrlpJunctIntronHsh{$junctStr}}) {
				my $hitSplicingType = ${$NGSJunctInfoHsh{$hitJunctStr}}{"splicingEffType"};
				if ((($splicingEffType eq "senseRare") and ($hitSplicingType eq "senseOften")) or (($splicingEffType eq "antisenseRare") and ($hitSplicingType eq "antisenseOften")) or (($hitSplicingType eq "senseRare") and ($splicingEffType eq "senseOften")) or (($hitSplicingType eq "antisenseRare") and ($splicingEffType eq "antisenseOften"))) {
					#----push the hit if the hit is in the same sense but opposite category
					push @tmpAry, $hitSplicingType;
				}
			}
		}
		
		if (@tmpAry == 0) {
			push @tmpAry, "null";
		} else {
			${$rareOftenOvrlpCountHsh{$superCtgry}}{$splicingEffType}++;
		}
		
		${${$NGSJunctInfoHsh{$junctStr}}{"rareOftenOvrlp"}} = join ";", @tmpAry;
		
	}
	
	#----print the category Stats
	open (CTGRYJUNCTCOUNT, ">$outDir/junctInfo/categoryBasedJunctCount.txt");
	
	print CTGRYJUNCTCOUNT "\nGene Based Category Sense Antisense counts\n";
	print CTGRYJUNCTCOUNT "ctgry\ttotal\tsense\tantisense\n";
	foreach my $ctgry (sort {$a cmp $b} keys %ctgryHitCountHsh) {
		my $total = keys %{${$ctgryHitCountHsh{$ctgry}}{"total"}};
		my $sense = keys %{${$ctgryHitCountHsh{$ctgry}}{"sense"}};
		my $antisense = keys %{${$ctgryHitCountHsh{$ctgry}}{"antisense"}};
		my $superCtgry = $superCtgryDefineHsh{$ctgry};
		if (not defined $superCtgry) {
			print CTGRYJUNCTCOUNT "warning: $ctgry has an undefined superCtgry\n";
		}
		
		${$superCtgryGeneCountHsh{$superCtgry}}{"total"} += $total;
		${$superCtgryGeneCountHsh{$superCtgry}}{"sense"} += $sense;
		${$superCtgryGeneCountHsh{$superCtgry}}{"antisense"} += $antisense;
		
		print CTGRYJUNCTCOUNT $ctgry."\t".$total."\t".$sense."\t".$antisense."\n";
	}
	
	print CTGRYJUNCTCOUNT "\nGene Based SuperCategory Sense Antisense counts\n";
	print CTGRYJUNCTCOUNT "superCtgry\ttotal\tsense\tantisense\n";
	foreach my $superCtgry (sort {$a cmp $b} keys %superCtgryGeneCountHsh) {
		my $total = ${$superCtgryGeneCountHsh{$superCtgry}}{"total"};
		my $sense = ${$superCtgryGeneCountHsh{$superCtgry}}{"sense"};
		my $antisense = ${$superCtgryGeneCountHsh{$superCtgry}}{"antisense"};

		print CTGRYJUNCTCOUNT $superCtgry."\t".$total."\t".$sense."\t".$antisense."\n";
	}
	
	print CTGRYJUNCTCOUNT "\nGene Based Category Often Rare counts\n";
	print CTGRYJUNCTCOUNT "ctgry\ttotal\tRare\tOften\n";
	foreach my $ctgry (keys %ctgrySenseRegularSpliceGeneCountHsh) {
		my $total = keys %{${$ctgryHitCountHsh{$ctgry}}{"total"}};
		my $rareNum = keys %{${$ctgrySenseRegularSpliceGeneCountHsh{$ctgry}}{"rare"}};
		my $oftenNum = keys %{${$ctgrySenseRegularSpliceGeneCountHsh{$ctgry}}{"often"}};
		print CTGRYJUNCTCOUNT $ctgry."\t".$total."\t".$rareNum."\t".$oftenNum."\n";
	}

	print CTGRYJUNCTCOUNT "\nGene Based Super Category Regular Rare counts\n";
	print CTGRYJUNCTCOUNT "superCtgry\ttotal\tsenseRare\tsenseOften\tantisenseRare\tantisenseOften\n";
	foreach my $superCtgry (keys %superCtgrySenseRegularSpliceGeneCountHsh) {
		my $total = ${$superCtgryGeneCountHsh{$superCtgry}}{"total"};
		my $senseRare = keys %{${$superCtgrySenseRegularSpliceGeneCountHsh{$superCtgry}}{"senseRare"}};
		my $senseOften = keys %{${$superCtgrySenseRegularSpliceGeneCountHsh{$superCtgry}}{"senseOften"}};
		my $antisenseRare = keys %{${$superCtgrySenseRegularSpliceGeneCountHsh{$superCtgry}}{"antisenseRare"}};
		my $antisenseOften = keys %{${$superCtgrySenseRegularSpliceGeneCountHsh{$superCtgry}}{"antisenseOften"}};
		print CTGRYJUNCTCOUNT $superCtgry."\t".$total."\t".$senseRare."\t".$senseOften."\t".$antisenseRare."\t".$antisenseOften."\n";
	}
	
	print CTGRYJUNCTCOUNT "\nJunction based Category counts\n";
	foreach my $ctrgyCombination (sort {$ctgryCountHsh{$b} <=> $ctgryCountHsh{$a}} keys %ctgryCountHsh) {
		my $total = ${$ctgryCountHsh{$ctrgyCombination}}{"total"};
		my $sense = ${$ctgryCountHsh{$ctrgyCombination}}{"sense"};
		my $antisense = ${$ctgryCountHsh{$ctrgyCombination}}{"antisense"};
		my $null = ${$ctgryCountHsh{$ctrgyCombination}}{"null"};
		my $senseRare = ${$ctgryCountHsh{$ctrgyCombination}}{"senseRare"};
		my $senseOften = ${$ctgryCountHsh{$ctrgyCombination}}{"senseOften"};
		my $antisenseRare = ${$ctgryCountHsh{$ctrgyCombination}}{"antisenseRare"};
		my $antisenseOften = ${$ctgryCountHsh{$ctrgyCombination}}{"antisenseOften"};

		my $pct = sprintf "%.06f", 100*$total/$totalJunctNum;
		print CTGRYJUNCTCOUNT $ctrgyCombination."\t".$total."\t".$pct."\t".$sense."\t".$antisense."\t".$senseOften."\t".$senseRare."\t".$antisenseOften."\t".$antisenseRare."\t".$null."\n";
	}
	
	print CTGRYJUNCTCOUNT "\nJunction based superCategory counts\n";
	foreach my $superCtgry (sort {$superCtgryJunctCountHsh{$b} <=> $superCtgryJunctCountHsh{$a}} keys %superCtgryJunctCountHsh) {
		my $total = ${$superCtgryJunctCountHsh{$superCtgry}}{"total"};
		my $sense = ${$superCtgryJunctCountHsh{$superCtgry}}{"sense"};
		my $antisense = ${$superCtgryJunctCountHsh{$superCtgry}}{"antisense"};
		my $null = ${$superCtgryJunctCountHsh{$superCtgry}}{"null"};
		my $senseRare = ${$superCtgryJunctCountHsh{$superCtgry}}{"senseRare"};
		my $senseOften = ${$superCtgryJunctCountHsh{$superCtgry}}{"senseOften"};
		my $antisenseRare = ${$superCtgryJunctCountHsh{$superCtgry}}{"antisenseRare"};
		my $antisenseOften = ${$superCtgryJunctCountHsh{$superCtgry}}{"antisenseOften"};
		my $SROvlpSF = ${$rareOftenOvrlpCountHsh{$superCtgry}}{"senseRare"};
		my $SFOvlpSR = ${$rareOftenOvrlpCountHsh{$superCtgry}}{"senseOften"};
		my $AROvlpAF = ${$rareOftenOvrlpCountHsh{$superCtgry}}{"antisenseRare"};
		my $AFOvlpAR = ${$rareOftenOvrlpCountHsh{$superCtgry}}{"antisenseOften"};
		
		my $pct = sprintf "%.06f", 100*$total/$totalJunctNum;
		print CTGRYJUNCTCOUNT $superCtgry."\t".$total."\t".$pct."\t".$sense."\t".$antisense."\t".$senseOften.":".$SFOvlpSR."\t".$senseRare.":".$SROvlpSF."\t".$antisenseOften.":".$AFOvlpAR."\t".$antisenseRare.":".$AROvlpAF."\t".$null."\n";
	}
	
	return (\%NGSJunctInfoHsh, $oftenSplicingEffCutoff, \%geneBasedJunctCountHsh);
	
}
########################################################################## calculateStandardDeviationAndMean
sub calculateStandardDeviationAndMean {

	my @numbers = @{$_[0]};

	my ($mean1, $std_dev);

	#Prevent division by 0 error in case you get junk data
	return undef unless(scalar(@numbers));

	# Step 1, find the mean of the numbers
	my $total1 = 0;
	foreach my $num (@numbers) {
		$total1 += $num;
	}
	$mean1 = $total1 / (scalar @numbers);
	
	if (@numbers >= 3) {
		
		# Step 2, find the mean of the squares of the differences
		# between each number and the mean
		my $total2 = 0;
		foreach my $num (@numbers) {
			$total2 += ($mean1-$num)**2;
		}
		my $mean2 = $total2 / (scalar @numbers);
	
		# Step 3, standard deviation is the square root of the
		# above mean
		$std_dev = sqrt($mean2);
	
	} else {
		
		$std_dev = "numberLessThan5"
	}

	return ($mean1, $std_dev);
}
########################################################################## printArrayContent
sub printArrayContent {
	
	my @aryToPrintAry = @{$_[0]};
	my $outputPath = $_[1];
	
	print "Printing $outputPath\n";
	
	open (ARYTOPRINTFILE, ">$outputPath");
	foreach my $lineToPrint (@aryToPrintAry) {
		print ARYTOPRINTFILE $lineToPrint."\n";
	}
	close ARYTOPRINTFILE;
	
}
########################################################################## statisticsOfMultiIntronTranscripts
sub statisticsOfMultiIntronTranscripts {

	my %refGeneIDByJunctStrHsh = %{$_[0]};
	my %altTrnscptInfoHsh = %{$_[1]};
	my %altJunctStrInfoHsh = %{$_[2]};
	my %junctReadNumHsh = %{$_[3]};
	my %countFturCovInfoHsh = %{$_[4]};
	my %NGSJunctInfoHsh = %{$_[5]};
	my %jClusterInfoHsh = %{$_[6]};
	my $minPrmntReadNumRatioAsNonStochasticAltSite = $_[7];
	
	print "Calculating the statistics on the intron Sizes.\n";
	
	#---get the distance of all possible combinations
	my %tmpJunstrByGeneIDHsh;
	foreach my $junctStr (keys %refGeneIDByJunctStrHsh) {
		my $geneID = $refGeneIDByJunctStrHsh{$junctStr};
		push @{$tmpJunstrByGeneIDHsh{$geneID}}, $junctStr;
	}
	
	my (@allPossibleSuperIntronSizeAry, @actualSuperIntronSizeAry, @allInferiorIntronSizeAry, @allRefJunctReadNumAry, @skippedJunctReadNumAry, @nonSkippedJunctReadNumAry, %superSkippedReadNumHsh, @skippedGeneCovNtAry, @nonSkippedGeneCovNtAry);
	my %exactSkippedRefJunctStrHsh;
	open (SUPERSKIPReadNum, ">$outDir/exonSkip/superSkippedReadNum.log.txt");
	open (NONSKIPPOSSIBLESUPERJReadNum, ">$outDir/exonSkip/nonSkippedPossibleSuperJunctReadNum.log.txt");
	open (MULTIINTRON, ">$outDir/splicingEfficienyOfAltSiteInMultiIntronTranscript.log.txt");
	print MULTIINTRON join "", (join "\t", ("geneID", "geneIntronNum", "refJunct", "junctStr", "readNum", "score", "strnd", "unq", "majorJunct", "prmntInCluster", "onPrmntIsofm", "superJunct", "hitRef", "clstJunctNum", "prmntReadNumRatio", "majorCluster", "splicingRatio", "geneCovNt", "senseSplicingEfficiency", "splicingEffType", "prmntDiff5", "prmntDiff3", "prmntDiffBoth", "intronSize", "intronFrame", "refOrAlt", "isofmID", "altIsofmFullORF", "altSplcType", "onRefGene", "refJunctAvgReadNum", "NGSRefJunctReadNumRatio", "readNumCovPerNtRatio", "altSiteType", "ovrlpRefJunctRdRatio", "avgOvrlpRefJReadNum")), "\n";

	my %nonStochasticAltSiteHsh;

	foreach my $geneID (keys %tmpJunstrByGeneIDHsh) {

		my $geneCovNt = ${$countFturCovInfoHsh{$geneID}}{"bothCov"}/${$countFturCovInfoHsh{$geneID}}{"length"};
		my $exonSkipped = "no";
		
		if (@{$tmpJunstrByGeneIDHsh{$geneID}} > 1) {#---genes with more than 1 intron
			my $geneIntronNum = @{$tmpJunstrByGeneIDHsh{$geneID}};
			my @tmpJunctBoundAry;
			my $cntg;
			foreach my $junctStr (@{$tmpJunstrByGeneIDHsh{$geneID}}) {
				my $refJunct = $junctStr;
				my $readNum = $junctReadNumHsh{$junctStr};
				$readNum = 0 if (not defined $readNum);
				push @allRefJunctReadNumAry, $readNum;
				my @junctStrSplt = split /:/, $junctStr;
				$cntg = $junctStrSplt[0];
				my $intronSize = $junctStrSplt[2] - $junctStrSplt[1] + 1;
				push @allInferiorIntronSizeAry, $intronSize;
				push @tmpJunctBoundAry, ($junctStrSplt[1], $junctStrSplt[2]);
				my $nonStochasticAltSite = 'no';
				
				#-----check for alternative site
				if (exists $NGSJunctInfoHsh{$junctStr}) {#---NGS confirmed
					my $cluster = ${$NGSJunctInfoHsh{$junctStr}}{'cluster'};
					foreach my $altJunctStr (@{${$jClusterInfoHsh{$cluster}}{"allJunctStr"}}) {
						my $readNum = ${$NGSJunctInfoHsh{$altJunctStr}}{"readNum"};
						my $score = ${$NGSJunctInfoHsh{$altJunctStr}}{"score"};
						my $strnd = ${$NGSJunctInfoHsh{$altJunctStr}}{"strnd"};
						my $unq = ${$NGSJunctInfoHsh{$altJunctStr}}{"unq"};
						my $majorJunct = ${$NGSJunctInfoHsh{$altJunctStr}}{"majorJunct"};
						my $prmntInCluster = ${$NGSJunctInfoHsh{$altJunctStr}}{"prmntInCluster"};
						my $onPrmntIsofm = ${$NGSJunctInfoHsh{$altJunctStr}}{"onPrmntIsofm"};
						my $superJunct = ${$NGSJunctInfoHsh{$altJunctStr}}{"superJunct"};
						my $hitRef = ${$NGSJunctInfoHsh{$altJunctStr}}{"hitRef"};
						my $clstJunctNum = ${$NGSJunctInfoHsh{$altJunctStr}}{"clstJunctNum"};
						my $prmntReadNumRatio = ${$NGSJunctInfoHsh{$altJunctStr}}{"prmntReadNumRatio"};
						my $majorCluster = ${$NGSJunctInfoHsh{$altJunctStr}}{"majorCluster"};
						my $splicingRatio = ${$NGSJunctInfoHsh{$altJunctStr}}{"splicingRatio"};
						my $geneCovNt = ${$NGSJunctInfoHsh{$altJunctStr}}{"geneCovNt"};
						my $senseSplicingEfficiency = ${$NGSJunctInfoHsh{$altJunctStr}}{"senseSplicingEfficiency"};
						my $splicingEffType = ${$NGSJunctInfoHsh{$altJunctStr}}{"splicingEffType"};
						my $prmntDiff5 = ${$NGSJunctInfoHsh{$altJunctStr}}{"prmntDiff5"};
						my $prmntDiff3 = ${$NGSJunctInfoHsh{$altJunctStr}}{"prmntDiff3"};
						my $prmntDiffBoth = ${$NGSJunctInfoHsh{$altJunctStr}}{"prmntDiffBoth"};
						my $intronSize = ${$NGSJunctInfoHsh{$altJunctStr}}{"intronSize"};
						my $intronFrame = ${$NGSJunctInfoHsh{$altJunctStr}}{"intronFrame"};
						my $refOrAlt = ${$NGSJunctInfoHsh{$altJunctStr}}{"refOrAlt"};
						my $isofmID = ${$NGSJunctInfoHsh{$altJunctStr}}{"isofmID"};
						my $altIsofmFullORF = ${$NGSJunctInfoHsh{$altJunctStr}}{"altIsofmFullORF"};
						my $altSplcType = ${$NGSJunctInfoHsh{$altJunctStr}}{"altSplcType"};
						my $onRefGene = ${$NGSJunctInfoHsh{$altJunctStr}}{"onRefGene"};
						my $refJunctAvgReadNum = ${$NGSJunctInfoHsh{$altJunctStr}}{"refJunctAvgReadNum"};
						my $NGSRefJunctReadNumRatio = ${$NGSJunctInfoHsh{$altJunctStr}}{"NGSRefJunctReadNumRatio"};
						my $readNumCovPerNtRatio = ${$NGSJunctInfoHsh{$altJunctStr}}{"readNumCovPerNtRatio"};
						my $altSiteType = ${$NGSJunctInfoHsh{$altJunctStr}}{"altSiteType"};
						my $ovrlpRefJunctRdRatio = ${$NGSJunctInfoHsh{$altJunctStr}}{"ovrlpRefJunctRdRatio"};
						my $avgOvrlpRefJReadNum = ${$NGSJunctInfoHsh{$altJunctStr}}{"avgOvrlpRefJReadNum"};
						
						if ($refOrAlt eq 'alt') {
							$nonStochasticAltSite = 'yes' if $prmntReadNumRatio >= $minPrmntReadNumRatioAsNonStochasticAltSite;
						}
						print MULTIINTRON join "", (join "\t", ($geneID, $geneIntronNum, $refJunct, $altJunctStr, $readNum, $score, $strnd, $unq, $majorJunct, $prmntInCluster, $onPrmntIsofm, $superJunct, $hitRef, $clstJunctNum, $prmntReadNumRatio, $majorCluster, $splicingRatio, $geneCovNt, $senseSplicingEfficiency, $splicingEffType, $prmntDiff5, $prmntDiff3, $prmntDiffBoth, $intronSize, $intronFrame, $refOrAlt, $isofmID, $altIsofmFullORF, $altSplcType, $onRefGene, $refJunctAvgReadNum, $NGSRefJunctReadNumRatio, $readNumCovPerNtRatio, $altSiteType, $ovrlpRefJunctRdRatio, $avgOvrlpRefJReadNum)), "\n";
					}
				} else {
					print MULTIINTRON join "", (join "\t", ($geneID, $geneIntronNum, $refJunct, 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null', 'null')), "\n";
				}
				${$nonStochasticAltSiteHsh{$geneID}}{$junctStr} = $nonStochasticAltSite;
			}
			@tmpJunctBoundAry = sort {$a <=> $b} @tmpJunctBoundAry;
			
			for (my $i = 0; $i <= ($#tmpJunctBoundAry - 3); $i += 2) {#---0, 2, 4 if $iMax is 7
				for (my $j = $i + 3; $j <= $#tmpJunctBoundAry; $j += 2) {
					my $superJunctStart = $tmpJunctBoundAry[$i];
					my $superJunctEnd = $tmpJunctBoundAry[$j];
					my $superIntronSize = $superJunctEnd - $superJunctStart + 1;
					my $possibleSuperJunctStr = join ":", ($cntg, $superJunctStart, $superJunctEnd);
					push @allPossibleSuperIntronSizeAry, $superIntronSize;
					my @tmpSkippedJunctReadNum;
					foreach my $junctStr (@{$tmpJunstrByGeneIDHsh{$geneID}}) {
						my @junctStrSplt = split /:/, $junctStr;
						if (($junctStrSplt[1] eq $superJunctStart) or ($junctStrSplt[2] eq $superJunctEnd)) {
							if ((exists $altJunctStrInfoHsh{$possibleSuperJunctStr}) and (${$altJunctStrInfoHsh{$possibleSuperJunctStr}}{"exonSkipExactType"} eq "exact")) {
								$exonSkipped = "yes";
								$exactSkippedRefJunctStrHsh{$junctStr}++;
							}
							my $skippedJunctReadNum = $junctReadNumHsh{$junctStr};
							push @tmpSkippedJunctReadNum, $skippedJunctReadNum;
						}
					}
					
					#---in case if the refJunct themselved doesnt have supported reads
					$tmpSkippedJunctReadNum[0] = 0 if (not defined $tmpSkippedJunctReadNum[0]);
					$tmpSkippedJunctReadNum[1] = 0 if (not defined $tmpSkippedJunctReadNum[1]);

					my $meanSkippedJunctReadNum = ($tmpSkippedJunctReadNum[0]+$tmpSkippedJunctReadNum[1])/2;
					
					if (exists $altJunctStrInfoHsh{$possibleSuperJunctStr}) {
						if (${$altJunctStrInfoHsh{$possibleSuperJunctStr}}{"exonSkipExactType"} eq "exact") {
							push @actualSuperIntronSizeAry, $superIntronSize;
							my $superJunctReadNum = $junctReadNumHsh{$possibleSuperJunctStr};
							${$superSkippedReadNumHsh{$possibleSuperJunctStr}}{"superJunctReadNum"} = $superJunctReadNum;
							${$superSkippedReadNumHsh{$possibleSuperJunctStr}}{"meanSkippedJunctReadNum"} = $meanSkippedJunctReadNum;
							print SUPERSKIPReadNum $superJunctReadNum."\t".$meanSkippedJunctReadNum."\n";
						}
					} else {
						print NONSKIPPOSSIBLESUPERJReadNum $meanSkippedJunctReadNum."\n";
					}
				}
			}
			
			if ($exonSkipped eq "yes") {
				push @skippedGeneCovNtAry, $geneCovNt;
			} elsif ($exonSkipped eq "no") {
				push @nonSkippedGeneCovNtAry, $geneCovNt;
			}
			
			#---count the readNum of junctions of genes with > 1 junction 
			foreach my $junctStr (@{$tmpJunstrByGeneIDHsh{$geneID}}) {
				my $readNum = $junctReadNumHsh{$junctStr};
				$readNum = 0 if (not defined $readNum);
				push @allRefJunctReadNumAry, $readNum;
				if (exists $exactSkippedRefJunctStrHsh{$junctStr}) {
					push @skippedJunctReadNumAry, $readNum;
				} else {
					push @nonSkippedJunctReadNumAry, $readNum;
				}
			}
		}
	}
	close MULTIINTRON;
	open (NONSTOCHASTICALSITE, ">$outDir/nonStochasticAltSiteOnMultipleIntronGene.log.xls");
	print NONSTOCHASTICALSITE join '', (join "\t", ("geneID", "refJunctNum", "nonStochasticAltSiteNum")), "\n";
	foreach my $geneID (sort {$a cmp $b} keys %nonStochasticAltSiteHsh) {
		my $refJunctNum = keys %{$nonStochasticAltSiteHsh{$geneID}};
		my $nonStochasticAltSiteNum = 0;
		foreach my $junctStr (keys %{$nonStochasticAltSiteHsh{$geneID}}) {
			$nonStochasticAltSiteNum++ if ${$nonStochasticAltSiteHsh{$geneID}}{$junctStr} eq 'yes';
		}
		print NONSTOCHASTICALSITE join '', (join "\t", ($geneID, $refJunctNum, $nonStochasticAltSiteNum)), "\n";
	}
	close NONSTOCHASTICALSITE;
	
	printArrayContent(\@skippedGeneCovNtAry, "$outDir/exonSkip/exonSkippedGeneCovNt.log.txt");
	printArrayContent(\@nonSkippedGeneCovNtAry, "$outDir/exonSkip/nonExonSkippedGeneCovNt.log.txt");
	printArrayContent(\@allPossibleSuperIntronSizeAry, "$outDir/allPossibleSuperIntronSize.log.txt");
	printArrayContent(\@actualSuperIntronSizeAry, "$outDir/actualSuperIntronSize.log.txt");
	printArrayContent(\@allInferiorIntronSizeAry, "$outDir/allInferiorIntronSize.log.txt");
	printArrayContent(\@skippedJunctReadNumAry, "$outDir/exonSkip/skippedJunctReadNum.log.txt");
	printArrayContent(\@allRefJunctReadNumAry, "$outDir/allRefJunctReadNum.log.txt");
	printArrayContent(\@nonSkippedJunctReadNumAry, "$outDir/exonSkip/nonSkippedJunctReadNum.log.txt");
}
########################################################################## analyzeNonCanAndCanNGSSOverlapping
sub summarizeNonCanAndCanNGSSOverlapping {

	#---my ($nonCanNGSJunctInfoHsh_ref, $nonCanNGSJunctFilterStrndHsh_ref) = summarizeNonCanAndCanNGSSOverlapping($XSOvrlpNonCanNGSJCanNGSJHsh_ref, $nonCanNGSJunctStrndHsh_ref, $NGSJunctStrndHsh_ref, $nonCanNGSJunctCntgHsh_ref, $nonCanJunctScoreHsh_ref, $nonCanJunctReadNumHsh_ref);

	my %XSOvrlpNonCanNGSJCanNGSJHsh = %{$_[0]};
	my %nonCanNGSJunctDummyStrndHsh = %{$_[1]};
	my %NGSJunctStrndHsh = %{$_[2]};
	my %nonCanNGSJunctCntgHsh = %{$_[3]};
	my %nonCanJunctScoreHsh = %{$_[4]};
	my %nonCanJunctReadNumHsh = %{$_[5]};
	
	my %nonCanNGSJunctInfoHsh;
	my %nonCanNGSJunctStrndHsh;
	
	print "Summarizing overlapping of canonical and noncanoical junctions\n";
	
	my %ovrlpScenrioHsh;
	
	foreach my $nonCanNGSJunctStr (keys %nonCanNGSJunctDummyStrndHsh) {

		my $nonCanNGSJunctStrnd = "*";
		my $numOvrlpCanJunct = 0;
		my $readNum = $nonCanJunctReadNumHsh{$nonCanNGSJunctStr};
		my $score = $nonCanJunctScoreHsh{$nonCanNGSJunctStr};
		my $cntg = $nonCanNGSJunctCntgHsh{$nonCanNGSJunctStr};

		if (exists $XSOvrlpNonCanNGSJCanNGSJHsh{$nonCanNGSJunctStr}) {
			my %tmpStrndHsh;
			foreach my $canNGSJunctStr (keys %{$XSOvrlpNonCanNGSJCanNGSJHsh{$nonCanNGSJunctStr}}) {
				my $strand = $NGSJunctStrndHsh{$canNGSJunctStr};
				$nonCanNGSJunctStrnd = $strand; #--overwrite each time
				$tmpStrndHsh{$strand}++;
				$numOvrlpCanJunct++;
			}
			
			my $numStrand = keys %tmpStrndHsh;
			if ($numStrand == 1) {#---the nonCanonical junct is overlapping with only canJunct of 1 direction
				$nonCanNGSJunctStrndHsh{$nonCanNGSJunctStr} = $nonCanNGSJunctStrnd;
			} else {
				$nonCanNGSJunctStrnd = "*";
			}
		}
		
		#---record the scenerios
		my $canJOvrlpScenerio;
		if ($numOvrlpCanJunct < 1) {
			$canJOvrlpScenerio = "noCanJunct";
		} elsif ($numOvrlpCanJunct == 1) {
			$canJOvrlpScenerio = "singleCanJunctSingleStrnd";
		} elsif ($numOvrlpCanJunct > 1) {
			if (exists $nonCanNGSJunctStrndHsh{$nonCanNGSJunctStr}) {
				$canJOvrlpScenerio = "multiCanJunctSingleStrnd";
			} else {
				$canJOvrlpScenerio = "multiCanJunctBothStrnd";
			}
		} else {
			die "debug: impossible scenerio\n";
		}
		$ovrlpScenrioHsh{$canJOvrlpScenerio}++;

		${$nonCanNGSJunctInfoHsh{$nonCanNGSJunctStr}}{"canJOvrlpScenerio"} = $canJOvrlpScenerio;
		${$nonCanNGSJunctInfoHsh{$nonCanNGSJunctStr}}{"strnd"} = $nonCanNGSJunctStrnd;
		${$nonCanNGSJunctInfoHsh{$nonCanNGSJunctStr}}{"numOvrlpCanJunct"} = $numOvrlpCanJunct;
		${$nonCanNGSJunctInfoHsh{$nonCanNGSJunctStr}}{"readNum"} = $readNum;
		${$nonCanNGSJunctInfoHsh{$nonCanNGSJunctStr}}{"score"} = $score;
		${$nonCanNGSJunctInfoHsh{$nonCanNGSJunctStr}}{"cntg"} = $cntg;
	}

	#---print the simple stats
	my $totalNonCanJunctNum = keys %nonCanNGSJunctDummyStrndHsh;
	foreach my $canJOvrlpScenerio (sort {$ovrlpScenrioHsh{$b} <=> $ovrlpScenrioHsh{$a}} keys %ovrlpScenrioHsh) {
		my $count = $ovrlpScenrioHsh{$canJOvrlpScenerio};
		my $pct = sprintf "%.04f", 100*$count/$totalNonCanJunctNum;
		print $count."\t".$pct."\t".$canJOvrlpScenerio."\n";
	}
	
	return \%nonCanNGSJunctInfoHsh, \%nonCanNGSJunctStrndHsh;
	
}
########################################################################## filterAndPrintNonCanNGSJInfo
sub filterAndPrintNonCanNGSJInfo {

	#my $filternonCanNGSJunctFilterStrndHsh_ref = filterAndPrintNonCanNGSJInfo(\%nonCanNGSJunctInfoHsh, $maxUnq, $minScore, $minReadNum);
	
	my %nonCanNGSJunctInfoHsh = %{$_[0]};
	my $maxUnq = $_[1];
	my $minScore = $_[2];
	my $minReadNum = $_[3];
	my $maxLoCmPrptn = $_[4];
	
	print "Filtering non-canoical junctions information\n";
	
	open (NONCANNGSJINFO, ">$outDir/junctInfo/all.nonCanNGSJ.info.log.txt");
	print NONCANNGSJINFO join "", ((join "\t", ("junctStr", "cntg", "loCmPrptn", "unq", "strnd", "numOvrlpCanJunct", "readNum", "score", "canJOvrlpScenerio", "leftFlankSeq", "rightFlankSeq", "left2ntJnctnSeq", "right2ntJnctnSeq", "completeSeq")), "\n");
	open (FILTERNONCANNGSJINFO, ">$outDir/junctInfo/filter.nonCanNGSJ.info.log.txt");
	print FILTERNONCANNGSJINFO join "", ((join "\t", ("junctStr", "cntg", "loCmPrptn", "unq", "strnd", "numOvrlpCanJunct", "readNum", "score", "canJOvrlpScenerio", "leftFlankSeq", "rightFlankSeq", "left2ntJnctnSeq", "right2ntJnctnSeq", "completeSeq")), "\n");
	
	my %filterOvrlpScenrioHsh;
	my %filterNonCanNGSJunctStrndHsh;

	foreach my $junctStr (keys %nonCanNGSJunctInfoHsh) {

		my $canJOvrlpScenerio = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"canJOvrlpScenerio"};
		my $strnd = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"strnd"};
		my $numOvrlpCanJunct = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"numOvrlpCanJunct"};
		my $readNum = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"readNum"};
		my $score = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"score"};
		my $cntg = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"cntg"};
		my $loCmPrptn = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"loCmPrptn"};
		my $unq = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"unq"};
		my $leftFlankSeq = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"leftFlankSeq"};
		my $rightFlankSeq = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"rightFlankSeq"};
		my $left2ntJnctnSeq = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"upSite"};
		my $right2ntJnctnSeq = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"downSite"};
		my $completeSeq = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"completeSeq"};

		print NONCANNGSJINFO join "", ((join "\t", ($junctStr, $cntg, $loCmPrptn, $unq, $strnd, $numOvrlpCanJunct, $readNum, $score, $canJOvrlpScenerio, $leftFlankSeq, $rightFlankSeq, $left2ntJnctnSeq, $right2ntJnctnSeq, $completeSeq)), "\n");
		
		if (($unq <= $maxUnq) and ($score >= $minScore) and ($readNum >= $minReadNum) and ($loCmPrptn <= $maxLoCmPrptn)){
			$filterOvrlpScenrioHsh{$canJOvrlpScenerio}++;
			$filterNonCanNGSJunctStrndHsh{$junctStr} = $strnd;
			print FILTERNONCANNGSJINFO join "", ((join "\t", ($junctStr, $cntg, $loCmPrptn, $unq, $strnd, $numOvrlpCanJunct, $readNum, $score, $canJOvrlpScenerio, $leftFlankSeq, $rightFlankSeq, $left2ntJnctnSeq, $right2ntJnctnSeq, $completeSeq)), "\n");
		}
	}
	
	
	my $totalFilterNonCanJunctNum = keys %filterNonCanNGSJunctStrndHsh;
	foreach my $canJOvrlpScenerio (sort {$filterOvrlpScenrioHsh{$b} <=> $filterOvrlpScenrioHsh{$a}} keys %filterOvrlpScenrioHsh) {
		my $count = $filterOvrlpScenrioHsh{$canJOvrlpScenerio};
		my $pct = sprintf "%.04f", 100*$count/$totalFilterNonCanJunctNum;
		print $count."\t".$pct."\t".$canJOvrlpScenerio."\n";
	}

	close NONCANNGSJINFO;
	
	return \%filterNonCanNGSJunctStrndHsh;
}
########################################################################## summarizeBothCanNonCanNGSJunctAltSiteOnRefJunct
sub summarizeBothCanNonCanNGSJunctAltSiteOnRefJunct {

	#---summarizeBothCanNonCanNGSJunctAltSiteOnRefJunct($SSOvrlpRefJNGSJHsh_ref, $SSOvrlpRefJNonCanNGSJHsh_ref, \%NGSJunctInfoHsh, \%refJunctInfoHsh, \%nonCanNGSJunctInfoHsh, $maxUnq, $maxLoCmPrptn);

	my %SSOvrlpNonCanNGSJRefTrnscptHsh = %{$_[0]};
	my %SSOvrlpRefJNGSJHsh = %{$_[1]};
	my %SSOvrlpRefJNonCanNGSJHsh = %{$_[2]};
	my %NGSJunctInfoHsh = %{$_[3]};
	my %refJunctInfoHsh = %{$_[4]};
	my %nonCanNGSJunctInfoHsh = %{$_[5]};
	my $maxUnq = $_[6];
	my $maxLoCmPrptn = $_[7];
	my $maxPolyBasePrptn = $_[8];
	
	my %canRefAltSiteInfoHsh;
	my %nonCanRefAltSiteInfoHsh;
	my %cnfrmdRefJunctInfoHsh;
	
	print "Summarizing the canonical and non-canonical alternative junctions on confirmed reference junctions\n";
	
	foreach my $refJunctStr  (keys %refJunctInfoHsh) {
		
		next if (not exists $NGSJunctInfoHsh{$refJunctStr}); #----the refJ is not confirmed by NGSJ
		
		my $loCmPrptn = ${$refJunctInfoHsh{$refJunctStr}}{"loCmPrptn"};
		my $unq = ${$refJunctInfoHsh{$refJunctStr}}{"unq"};
		my $refLeftFlankSeq = ${$refJunctInfoHsh{$refJunctStr}}{"leftFlankSeq"};
		my $refRightFlankSeq = ${$refJunctInfoHsh{$refJunctStr}}{"rightFlankSeq"};
		my $refLeft2ntJnctnSeq = ${$refJunctInfoHsh{$refJunctStr}}{"upSite"};
		my $refRight2ntJnctnSeq = ${$refJunctInfoHsh{$refJunctStr}}{"downSite"};
		my $strnd = ${$refJunctInfoHsh{$refJunctStr}}{"strnd"};

		my @refJunctStrSplt = split /:/, $refJunctStr;
		my $refJunctBoundStart = $refJunctStrSplt[1]-1;
		my $refJunctBoundEnd = $refJunctStrSplt[2]+1;
		
		#----skip the junctions in low complexity region and 
		if (($unq <= $maxUnq) and ($loCmPrptn <= $maxLoCmPrptn)) {

			#---overlapping canonical junctions
			if (exists $SSOvrlpRefJNGSJHsh{$refJunctStr}) {
				foreach my $junctStr (keys %{$SSOvrlpRefJNGSJHsh{$refJunctStr}}) {
					my $readNum = ${$NGSJunctInfoHsh{$junctStr}}{"readNum"};
					my $score = ${$NGSJunctInfoHsh{$junctStr}}{"score"};
					my $leftFlankSeq = ${$NGSJunctInfoHsh{$junctStr}}{"leftFlankSeq"};
					my $rightFlankSeq = ${$NGSJunctInfoHsh{$junctStr}}{"rightFlankSeq"};
					my $upSite = ${$NGSJunctInfoHsh{$junctStr}}{"upSite"};
					my $downSite = ${$NGSJunctInfoHsh{$junctStr}}{"downSite"};
					my $superJunct = ${$NGSJunctInfoHsh{$junctStr}}{"superJunct"};
					my $splicingRatio = ${$NGSJunctInfoHsh{$junctStr}}{"splicingRatio"};
					my $senseSplicingEfficiency = ${$NGSJunctInfoHsh{$junctStr}}{"senseSplicingEfficiency"};
					my $splicingEffType = ${$NGSJunctInfoHsh{$junctStr}}{"splicingEffType"};
					my $refOrAlt = ${$NGSJunctInfoHsh{$junctStr}}{"refOrAlt"};
					my $leftExonPolyBasePrptn = ${$NGSJunctInfoHsh{$junctStr}}{"leftExonPolyBasePrptn"};
					my $rightExonPolyBasePrptn = ${$NGSJunctInfoHsh{$junctStr}}{"rightExonPolyBasePrptn"};
					my $leftExonSeq = ${$NGSJunctInfoHsh{$junctStr}}{"leftExonSeq"};
					my $rightExonSeq = ${$NGSJunctInfoHsh{$junctStr}}{"rightExonSeq"};
					
					my @junctStrSplt = split /:/, $junctStr;
					my $NGSJunctBoundStart = $junctStrSplt[1]-1;
					my $NGSJunctBoundEnd = $junctStrSplt[2]+1;
					
					my $refShift5 = "refShift5";
					my $refShift3 = "refShift3";
					
					if ($strnd eq "+") {
						$refShift5 = $NGSJunctBoundStart - $refJunctBoundStart;
						$refShift3 = $NGSJunctBoundEnd - $refJunctBoundEnd;
					} elsif ($strnd eq "-") {
						$refShift3 = -1*($NGSJunctBoundStart - $refJunctBoundStart);
						$refShift5 = -1*($NGSJunctBoundEnd - $refJunctBoundEnd);
					} else {
						die "impossible strand\n";
					}

					if ($refOrAlt eq "ref") {#---this is the reference junction

						${$cnfrmdRefJunctInfoHsh{$refJunctStr}}{"readNum"} = $readNum;
						${$cnfrmdRefJunctInfoHsh{$refJunctStr}}{"score"} = $score;
						${$cnfrmdRefJunctInfoHsh{$refJunctStr}}{"strnd"} = $strnd;
						${$cnfrmdRefJunctInfoHsh{$refJunctStr}}{"splicingEffType"} = $splicingEffType;
						${$cnfrmdRefJunctInfoHsh{$refJunctStr}}{"senseSplicingEfficiency"} = $senseSplicingEfficiency;
						${$cnfrmdRefJunctInfoHsh{$refJunctStr}}{"splicingRatio"} = $splicingRatio;
						
					} elsif ($refOrAlt eq "alt") {#---alt

						if ($superJunct ne "y") {#---dont count super junctions
							if (($leftExonPolyBasePrptn <= $maxPolyBasePrptn) and ($rightExonPolyBasePrptn <= $maxPolyBasePrptn) and ($refShift5 != $refShift3)) {

								${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"readNum"} = $readNum;
								${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"score"} = $score;
								${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"splicingEffType"} = $splicingEffType;
								${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"senseSplicingEfficiency"} = $senseSplicingEfficiency;
								${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"splicingRatio"} = $splicingRatio;
								${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"leftFlankSeq"} = $leftFlankSeq;
								${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"rightFlankSeq"} = $rightFlankSeq;
								${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"upSite"} = $upSite;
								${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"downSite"} = $downSite;
								${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"leftExonPolyBasePrptn"} = $leftExonPolyBasePrptn;
								${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"rightExonPolyBasePrptn"} = $rightExonPolyBasePrptn;
								${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"leftExonSeq"} = $leftExonSeq;
								${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"rightExonSeq"} = $rightExonSeq;
								${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"refShift5"} = $refShift5;
								${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"refShift3"} = $refShift3;
							}
						}

					} else {#--others, like long juncts cross multiple genes
						#--do nothing at the moment
					}
				}
			}
			
			#---overlapping non-canonical junctions
			if (exists $SSOvrlpRefJNonCanNGSJHsh{$refJunctStr}) {
				foreach my $junctStr (keys %{$SSOvrlpRefJNonCanNGSJHsh{$refJunctStr}}) {

					if ((keys %{$SSOvrlpNonCanNGSJRefTrnscptHsh{$junctStr}}) == 1) {#---overlap with only one refTranscript
						
						my $ovrlpScenerio;
						foreach my $refTrnscpt (keys %{$SSOvrlpNonCanNGSJRefTrnscptHsh{$junctStr}}) {
							$ovrlpScenerio = ${$SSOvrlpNonCanNGSJRefTrnscptHsh{$junctStr}}{$refTrnscpt};
						}

						if ($ovrlpScenerio == 3) {#--- junction is within the gene

							my @junctStrSplt = split /:/, $junctStr;
							my $NGSJunctBoundStart = $junctStrSplt[1]-1;
							my $NGSJunctBoundEnd = $junctStrSplt[2]+1;
							
							my $refShift5 = "refShift5";
							my $refShift3 = "refShift3";
							
							if ($strnd eq "+") {
								$refShift5 = $NGSJunctBoundStart - $refJunctBoundStart;
								$refShift3 = $NGSJunctBoundEnd - $refJunctBoundEnd;
							} elsif ($strnd eq "-") {
								$refShift3 = -1*($NGSJunctBoundStart - $refJunctBoundStart);
								$refShift5 = -1*($NGSJunctBoundEnd - $refJunctBoundEnd);
							} else {
								die "impossible strand\n";
							}

							my $numOvrlpCanJunct = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"numOvrlpCanJunct"};
							my $readNum = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"readNum"};
							my $score = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"score"};
							my $leftFlankSeq = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"leftFlankSeq"};
							my $rightFlankSeq = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"rightFlankSeq"};
							my $upSite = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"upSite"};
							my $downSite = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"downSite"};
							my $senseSplicingEfficiency = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"senseSplicingEfficiency"};
							my $splicingEffType = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"splicingEffType"};
							my $leftExonPolyBasePrptn = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"leftExonPolyBasePrptn"};
							my $rightExonPolyBasePrptn = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"rightExonPolyBasePrptn"};
							my $leftExonSeq = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"leftExonSeq"};
							my $rightExonSeq = ${$nonCanNGSJunctInfoHsh{$junctStr}}{"rightExonSeq"};
							
							if (($leftExonPolyBasePrptn <= $maxPolyBasePrptn) and ($rightExonPolyBasePrptn <= $maxPolyBasePrptn) and ($refShift5 != $refShift3)) {
								${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"readNum"} = $readNum;
								${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"score"} = $score;
								${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"splicingEffType"} = $splicingEffType;
								${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"senseSplicingEfficiency"} = $senseSplicingEfficiency;
								${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"leftFlankSeq"} = $leftFlankSeq;
								${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"rightFlankSeq"} = $rightFlankSeq;
								${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"upSite"} = $upSite;
								${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"downSite"} = $downSite;
								${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"leftExonPolyBasePrptn"} = $leftExonPolyBasePrptn;
								${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"rightExonPolyBasePrptn"} = $rightExonPolyBasePrptn;
								${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"leftExonSeq"} = $leftExonSeq;
								${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"rightExonSeq"} = $rightExonSeq;
								${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"refShift5"} = $refShift5;
								${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"refShift3"} = $refShift3;
							}
						}#---if ($ovrlpScenerio == 3) {#--- junction is within the gene	
					} #---if ((keys $SSOvrlpNonCanNGSJRefTrnscptHsh{$junctStr}) == 1) {#---overlap with only one refTranscript
				}#---foreach my $junctStr (keys %{$SSOvrlpRefJNonCanNGSJHsh{$refJunctStr}})
			}#---if (exists $SSOvrlpRefJNonCanNGSJHsh{$refJunctStr})
		}#---if (($unq <= $maxUnq) and ($loCmPrptn <= $maxLoCmPrptn))
	}#---foreach my $refJunctStr  (keys %refJunctInfoHsh)
	
	open (REFJALTINFO, ">$outDir/junctInfo/refJunct.altSite.info.log.txt");
	open (ALTJALTINFO, ">$outDir/junctInfo/altJunct.altSite.info.log.txt");
	print ALTJALTINFO join "", ((join "\t", ("canonicality", "refJunctStr", "junctStr", "score", "readNum", "refRdNumRatio", "splicingEffType", "senseSplicingEfficiency", "leftFlankSeq", "rightFlankSeq", "upSite", "downSite", "refShift3", "refShift5", "leftExonPolyBasePrptn", "rightExonPolyBasePrptn", "leftExonSeq", "rightExonSeq")), "\n");
	
	foreach my $refJunctStr (keys %cnfrmdRefJunctInfoHsh) {
		
		my $refReadNum = ${$cnfrmdRefJunctInfoHsh{$refJunctStr}}{"readNum"};
		my $refScore = ${$cnfrmdRefJunctInfoHsh{$refJunctStr}}{"score"};
		my $strnd = ${$cnfrmdRefJunctInfoHsh{$refJunctStr}}{"strnd"};
		my $refSplicingEffType = ${$cnfrmdRefJunctInfoHsh{$refJunctStr}}{"splicingEffType"};
		my $refSenseSplicingEfficiency = ${$cnfrmdRefJunctInfoHsh{$refJunctStr}}{"senseSplicingEfficiency"};
		
		my $canAltJNum = my $nonCanAltJNum = my $canAltRdNum = my $nonCanAltRdNum = 0;

		if (exists $canRefAltSiteInfoHsh{$refJunctStr}) {
			foreach my $junctStr (keys %{$canRefAltSiteInfoHsh{$refJunctStr}}) {

				${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"refRdNumRatio"} = sprintf "%.06f", ${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"readNum"}/$refReadNum;

				my $readNum = ${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"readNum"};
				my $refRdNumRatio = ${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"refRdNumRatio"};
				my $score = ${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"score"};
				my $splicingEffType = ${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"splicingEffType"};
				my $senseSplicingEfficiency = ${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"senseSplicingEfficiency"};
				my $leftFlankSeq = ${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"leftFlankSeq"};
				my $rightFlankSeq = ${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"rightFlankSeq"};
				my $upSite = ${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"upSite"};
				my $downSite = ${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"downSite"};
				my $refShift3 = ${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"refShift3"};
				my $refShift5 = ${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"refShift5"};
				my $leftExonPolyBasePrptn = ${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"leftExonPolyBasePrptn"};
				my $rightExonPolyBasePrptn = ${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"rightExonPolyBasePrptn"};
				my $leftExonSeq = ${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"leftExonSeq"};
				my $rightExonSeq = ${${$canRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"rightExonSeq"};

				$canAltJNum++;
				$canAltRdNum += $readNum;
				
				print ALTJALTINFO join "", ((join "\t", ("c", $refJunctStr, $junctStr, $score, $readNum, $refRdNumRatio, $splicingEffType, $senseSplicingEfficiency, $leftFlankSeq, $rightFlankSeq, $upSite, $downSite, $refShift3, $refShift5, $leftExonPolyBasePrptn, $rightExonPolyBasePrptn, $leftExonSeq, $rightExonSeq)), "\n");
				
			}
		}
		
		if (exists $nonCanRefAltSiteInfoHsh{$refJunctStr}) {
			foreach my $junctStr (keys %{$nonCanRefAltSiteInfoHsh{$refJunctStr}}) {

				${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"refRdNumRatio"} = sprintf "%.06f", ${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"readNum"}/$refReadNum;

				my $readNum = ${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"readNum"};
				my $refRdNumRatio = ${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"refRdNumRatio"};
				my $score = ${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"score"};
				my $splicingEffType = ${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"splicingEffType"};
				my $senseSplicingEfficiency = ${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"senseSplicingEfficiency"};
				my $leftFlankSeq = ${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"leftFlankSeq"};
				my $rightFlankSeq = ${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"rightFlankSeq"};
				my $upSite = ${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"upSite"};
				my $downSite = ${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"downSite"};
				my $refShift3 = ${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"refShift3"};
				my $refShift5 = ${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"refShift5"};
				my $leftExonPolyBasePrptn = ${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"leftExonPolyBasePrptn"};
				my $rightExonPolyBasePrptn = ${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"rightExonPolyBasePrptn"};
				my $leftExonSeq = ${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"leftExonSeq"};
				my $rightExonSeq = ${${$nonCanRefAltSiteInfoHsh{$refJunctStr}}{$junctStr}}{"rightExonSeq"};

				$nonCanAltJNum++;
				$nonCanAltRdNum += $readNum;
				
				print ALTJALTINFO join "", ((join "\t", ("n", $refJunctStr, $junctStr, $score, $readNum, $refRdNumRatio, $splicingEffType, $senseSplicingEfficiency, $leftFlankSeq, $rightFlankSeq, $upSite, $downSite, $refShift3, $refShift5, $leftExonPolyBasePrptn, $rightExonPolyBasePrptn, $leftExonSeq, $rightExonSeq)), "\n");
				
			}		
		}
	}
	
	close REFJALTINFO;
	close ALTJALTINFO;
}
########################################################################## countNGSJunctOverlapFtur
sub summarizeNonCanNGSJunctOverlapFtur {
	
	#---my $nonCannonCanNGSJunctInfoHsh_ref = summarizeNonCanNGSJunctOverlapFtur(\%nonCannonCanNGSJunctInfoHsh, $SSOvrlpNonCanNGSJFturCountHsh_ref, $XSOvrlpNonCanNGSJFturCountHsh_ref, $countFturCovInfoHsh_ref, $nonCanJunctReadNumHsh_ref, $oftenSplicingEffCutoff)
	
	my %nonCannonCanNGSJunctInfoHsh = %{$_[0]};
	my %SSOvrlpNonCanNGSJFturCountHsh = %{$_[1]};
	my %XSOvrlpNonCanNGSJFturCountHsh = %{$_[2]};
	my %countFturCovInfoHsh = %{$_[3]};
	my %nonCanJunctReadNumHsh = %{$_[4]};
	my $oftenSplicingEffCutoff = $_[5];
	
	print "Summarizing nonCanonical junctions on counting fturs\n";
	
	#---go through all juncts
	foreach my $junctStr (keys %nonCanNGSJunctInfoHsh) {
		my $readNum = $nonCanJunctReadNumHsh{$junctStr};
		${$nonCanNGSJunctInfoHsh{$junctStr}}{"senseSplicingEfficiency"} = "null";
		${$nonCanNGSJunctInfoHsh{$junctStr}}{"antisenseSplicingEfficiency"} = "null";
		${$nonCanNGSJunctInfoHsh{$junctStr}}{"splicingEffType"} = "null";
		my $splicingEfficiency = -1;
		my @ctrgyAry = my @senseAry = my @antisenseAry;
		my $splicingEffType = "null";

		if (exists 	$XSOvrlpNonCanNGSJFturCountHsh{$junctStr}) {#---overlap with something
			my $geneNum = keys %{$XSOvrlpNonCanNGSJFturCountHsh{$junctStr}};
			if ($geneNum == 1) {#---record junctions only overlap with one ftur
				foreach my $geneID (keys %{$XSOvrlpNonCanNGSJFturCountHsh{$junctStr}}) {
				
					my $geneCovNt = (${$countFturCovInfoHsh{$geneID}}{"plusCov"} + ${$countFturCovInfoHsh{$geneID}}{"minusCov"}) / ${$countFturCovInfoHsh{$geneID}}{"length"};
					my $splicingEfficiency = 999;
					$splicingEfficiency = $readNum/$geneCovNt if ($geneCovNt > 0);

					if (exists ${$SSOvrlpNonCanNGSJFturCountHsh{$junctStr}}{$geneID}) {
						${$nonCanNGSJunctInfoHsh{$junctStr}}{"senseSplicingEfficiency"} = $splicingEfficiency;
						if ($splicingEfficiency < $oftenSplicingEffCutoff) {
							$splicingEffType = "senseRare";
						} else {#----> greater than $oftenSplicingEffCutoff
							$splicingEffType = "senseOften";
						}
					} else {
						if ($splicingEfficiency < $oftenSplicingEffCutoff) {
							$splicingEffType = "antisenseRare";
						} else {#----> greater than $oftenSplicingEffCutoff
							$splicingEffType = "antisenseOften";
						}
						${$nonCanNGSJunctInfoHsh{$junctStr}}{"antisenseSplicingEfficiency"} = $splicingEfficiency;
					}

					${$nonCanNGSJunctInfoHsh{$junctStr}}{"splicingEffType"} = $splicingEffType;

				}#---foreach my $geneID (keys %{$XSOvrlpNonCanNGSJFturCountHsh{$junctStr}})
			}#---if ($geneNum == 1)
		}#---end of if (exists 	$XSOvrlpNonCanNGSJFturCountHsh{$junctStr}) {#---overlap with something
	}#---end of foreach my $junctStr (keys %nonCanNGSJunctInfoHsh)

	return \%nonCanNGSJunctInfoHsh;
	
}
########################################################################## countBaseContent
sub countBaseContent {

	my $seqToCount = $_[0];
	my @charToCount = @{$_[1]};
	
	my $seqLen = length ($seqToCount);
	my %prptnHsh = ();
	my %countHsh = ();
	my $maxCount = 0;
	my $maxPrptn = 0;
	my $maxChar = "null";
	
	foreach my $char (@charToCount) {
		my $tmpSeq = $seqToCount;
		my $count = ($tmpSeq =~ s/$char//g);
		$count = 0 if ($count eq "");
		my $prptn = sprintf "%.06f", $count/$seqLen;
		$prptnHsh{$char} = $prptn;
		$countHsh{$char} = $count;
		if ($count > $maxCount) {
			$maxCount = $count;
			$maxChar = $char;
			$maxPrptn = $prptn;
		}
	}
	
	return $maxChar, $maxPrptn, $maxCount, \%prptnHsh, \%countHsh;
	
}
########################################################################## plotSplicingEfficiencyVersusSiteEntropy
sub plotSplicingEfficiencyVersusSiteEntropy {
	
	#--- plotSplicingEfficiencyVersusSiteEntropy(\%NGSJunctInfoHsh, 100, $intronlessGeneCodingSeqHsh_ref);
	
	my %NGSJunctInfoHsh = %{$_[0]};
	my $numJunctInterval = $_[1];
	my %intronlessGeneCodingSeqHsh = %{$_[2]};

	my $minReadNum = 1;
	my $maxReadNum = 99999;
	my $length = 50;
	my $minSenseSplicingEfficiency = 0.0000000000000000001;
	my $maxSenseSplicingEfficiency = 1;
	my $binWidthLogSenseSplicingEfficiency = 0.2;
	
	my $trimExonEnd = 40;
	my $trimIntronEnd = 40;
	
	my %senseSplicingEfficiencySortingHsh;
	
	#---get the junctStr with with complete sequence, with unq == 1 and with >0 senseSplicingEfficiency
	foreach my $junctStr (keys %NGSJunctInfoHsh) {
		my $unq = ${$NGSJunctInfoHsh{$junctStr}}{"unq"};
		my $completeSeq = ${$NGSJunctInfoHsh{$junctStr}}{"completeSeq"};
		my $senseSplicingEfficiency = ${$NGSJunctInfoHsh{$junctStr}}{"senseSplicingEfficiency"};
		my $readNum = ${$NGSJunctInfoHsh{$junctStr}}{"readNum"};
		
		${$NGSJunctInfoHsh{$junctStr}}{"consensusValue"} = "null";
		if (($completeSeq eq "complete") and ($unq <= 1.8) and ($senseSplicingEfficiency ne "null") and ($readNum >= $minReadNum) and ($readNum <= $maxReadNum) and ($senseSplicingEfficiency >= $minSenseSplicingEfficiency)  and ($senseSplicingEfficiency <= $maxSenseSplicingEfficiency)) {
			$senseSplicingEfficiencySortingHsh{$junctStr} = $senseSplicingEfficiency;
		}
	}
	my $validJnctNum = keys %senseSplicingEfficiencySortingHsh;
	
	print "\n$validJnctNum are valid.\n";
	
	my %posFactorHsh;
	
	#for my $pos (1..10) {$posFactorHsh{$pos} = 1;} #---5' site
	for my $pos (13..17) {$posFactorHsh{$pos} = 1;} #---5' site
	for my $pos (18..22) {$posFactorHsh{$pos} = 1;} #---5' site
	for my $pos (23..31) {$posFactorHsh{$pos} = 1;} #---5' site
	for my $pos (32..32) {$posFactorHsh{$pos} = 1;} #---5' site
	#for my $pos (35..44) {$posFactorHsh{$pos} = 1;} #---5' site
	
	#for my $pos (1..8) {$posFactorHsh{$pos} = 1;} #---5' site
	#for my $pos (11..15) {$posFactorHsh{$pos} = 3;} #---5' site
	#for my $pos (16..18) {$posFactorHsh{$pos} = 2;} #---5' site
	#for my $pos (19..25) {$posFactorHsh{$pos} = 2;} #--3' site
	#for my $pos (26..26) {$posFactorHsh{$pos} = 3;} #--3' site
	#for my $pos (29..36) {$posFactorHsh{$pos} = 1;} #---3' site
	
	my $indivconsensusValueSplicingEfficiencyLogPath = "$outDir/splicingEffInterval/log/indiv.consensusValueSplicingEfficiency.log.txt";
	open (INDIVCNSCRLOG, ">$indivconsensusValueSplicingEfficiencyLogPath");

	my $intervalconsensusValueSplicingEfficiencyLogPath = "$outDir/splicingEffInterval/log/interval.consensusValueSplicingEfficiency.log.txt";
	open (INTRVCNSCRLOG, ">$intervalconsensusValueSplicingEfficiencyLogPath");

	my $binBySpEffconsensusValueSplicingEfficiencyLogPath = "$outDir/splicingEffInterval/log/binBySpEff.$binWidthLogSenseSplicingEfficiency.consensusValueSplicingEfficiency.log.txt";
	open (EFFBINCNSCRLOG, ">$binBySpEffconsensusValueSplicingEfficiencyLogPath");

	my $similaritySplicingEfficiencyLogPath = "$outDir/splicingEffInterval/log/similaritySplicingEfficiency.log.txt";
	open (SMLRTYLOG, ">$similaritySplicingEfficiencyLogPath");

	my $exonSeqFastaPath = "$outDir/dreme/fasta/oftenSplicedExonSeq.fasta";
	open (EXONSEQFASTA, ">$exonSeqFastaPath");

	my $seqForRefPSSMFastaPath = "$outDir/splicingEffInterval/fasta/seqForRefPSSM.fasta";
	open (REFPSSMSEQ, ">$seqForRefPSSMFastaPath");

	#---generate the ref PSSM
	print "Generating refPSSM\n";
	
=pod
	my $refPSSMJunctRemoveCount = 0;
	my $pct5Num = int (0.00*$validJnctNum);
	my $topPct5Num = $pct5Num;
	my $bottomPct5Num = $validJnctNum - $pct5Num;

	foreach my $junctStr (sort {$senseSplicingEfficiencySortingHsh{$b} <=> $senseSplicingEfficiencySortingHsh{$a}} keys %senseSplicingEfficiencySortingHsh) {
		$refPSSMJunctRemoveCount++;
		delete $senseSplicingEfficiencySortingHsh{$junctStr} if (($refPSSMJunctRemoveCount <= $topPct5Num) or ($refPSSMJunctRemoveCount >= $bottomPct5Num));
	}
	
	my $pct90ValidJnctNum = keys %senseSplicingEfficiencySortingHsh;
	print "\nTop and bottom $pct5Num junctions removed, $pct90ValidJnctNum jounction left.\n";
=cut

	my $topRefPSSMJunctPct = 5;
	my $topRefPSSMJunctNum = int (($topRefPSSMJunctPct/100)*$validJnctNum); #---top pct
	my $refPSSMJunctCounted = 0;
	
	my %refPSSMPosSeqHsh;
	
	foreach my $junctStr (sort {$senseSplicingEfficiencySortingHsh{$b} <=> $senseSplicingEfficiencySortingHsh{$a}} keys %senseSplicingEfficiencySortingHsh) {
		$refPSSMJunctCounted++;
		my $leftFlankSeq = ${$NGSJunctInfoHsh{$junctStr}}{"leftFlankSeq"};
		my $rightFlankSeq = ${$NGSJunctInfoHsh{$junctStr}}{"rightFlankSeq"};
		my $oriLen = length $leftFlankSeq;
		my $trimLen = $oriLen-$trimExonEnd-$trimIntronEnd;
		die "the trimmed the length of the flank sequence is smaller than zero" if ($trimLen <= 0) ;
		my $trimLeftFlankSeq = substr $leftFlankSeq, $trimExonEnd, $trimLen;
		my $trimRightFlankSeq = substr $rightFlankSeq, $trimIntronEnd, $trimLen;
		my $concateSeq = $trimLeftFlankSeq.$trimRightFlankSeq;
		my $senseSplicingEfficiency = ${$NGSJunctInfoHsh{$junctStr}}{"senseSplicingEfficiency"};
		my $splicingEffType = ${$NGSJunctInfoHsh{$junctStr}}{"splicingEffType"};
		my @concateSeqSplt = split //, $concateSeq;
		my $posSeq = "";
		foreach my $i (0..$#concateSeqSplt) {
			my $pos = $i + 1;
			$posSeq .= $concateSeqSplt[$i] if (exists $posFactorHsh{$pos});
		}
		$refPSSMPosSeqHsh{$junctStr} = $posSeq;
		print REFPSSMSEQ ">$junctStr\n";
		print REFPSSMSEQ "$concateSeq\n";
		last if $refPSSMJunctCounted >= $topRefPSSMJunctNum;
	}
	close REFPSSMSEQ;
	
	my $refPSSMEntropyPath = "$outDir/splicingEffInterval/entropy/refPSSM.entropy.txt";
	my $refPSSMLogoPath = "$outDir/splicingEffInterval/pdf/refPSSM.pdf";
	
	my $refPSSMNumSeq = `sed -n \'\$=\' $seqForRefPSSMFastaPath`;
	$refPSSMNumSeq = $refPSSMNumSeq/2;
	my $refPSSMTitle = "oftenSpliced";
	$refPSSMTitle .= "_n=$refPSSMNumSeq";
	my $refPSSMWeblogoEntropyCMD = "weblogo -f $seqForRefPSSMFastaPath -F logodata -n $length -s large --title $refPSSMTitle -A rna -c classic >$refPSSMEntropyPath";
	my $refPSSMGrepCMD = "ps -ef | grep weblogo | grep -v grep | grep -v perl";
	runAndCheckSerialTask($refPSSMGrepCMD, "weblogo", $refPSSMWeblogoEntropyCMD, "$outDir/error.log.txt");

	$refPSSMWeblogoEntropyCMD = "weblogo -f $seqForRefPSSMFastaPath -F pdf -n $length -s large --title $refPSSMTitle -A rna -c classic >$refPSSMLogoPath";
	$refPSSMGrepCMD = "ps -ef | grep weblogo | grep -v grep | grep -v perl";
	runAndCheckSerialTask($refPSSMGrepCMD, "weblogo", $refPSSMWeblogoEntropyCMD, "$outDir/error.log.txt");
	
	#----read the refPSSM
	my %refPSSMHsh;
	open (REFENTROPY, "$refPSSMEntropyPath");
	while (my $theLine = <REFENTROPY>) {
		chomp $theLine;
		$theLine =~ s/ //g;
		next if (($theLine =~ m/^#/) or (length $theLine < 2));
		my ($pos, $ACount, $CCount, $GCount, $UCount, $meanEntropy, $lowEntropy, $highEntropy, $weight) = split /\t/, $theLine;
		if (exists $posFactorHsh{$pos}) {
			my $factor = $posFactorHsh{$pos};
			my $totalCount = $ACount + $CCount + $GCount + $UCount;
			#my $totalCount = 1; #----switch off proptn
			my @sortCountAry = sort {$a <=> $b} ($ACount, $CCount, $GCount, $UCount);
			${$refPSSMHsh{$pos}}{"max"} = $sortCountAry[-1]/$totalCount;
			${$refPSSMHsh{$pos}}{"min"} = $sortCountAry[0]/$totalCount;
			${$refPSSMHsh{$pos}}{"A"} = $ACount/$totalCount;
			${$refPSSMHsh{$pos}}{"C"} = $CCount/$totalCount;
			${$refPSSMHsh{$pos}}{"G"} = $GCount/$totalCount;
			${$refPSSMHsh{$pos}}{"U"} = $UCount/$totalCount;
			${$refPSSMHsh{$pos}}{"meanEntropy"} = $meanEntropy;
			${$refPSSMHsh{$pos}}{"lowEntropy"} = $lowEntropy;
			${$refPSSMHsh{$pos}}{"highEntropy"} = $highEntropy;
			${$refPSSMHsh{$pos}}{"factor"} = $factor;
		}
	}
	close REFENTROPY;

	open (REFENTROPY, "$refPSSMEntropyPath");
	
	my $upperbinLimit = log($maxSenseSplicingEfficiency)/log(10);
	my $lowerbinLimit = (log($maxSenseSplicingEfficiency)/log(10)) - $binWidthLogSenseSplicingEfficiency;
	my (@tmpJunctIntrvlSplcngEffAry, @tmpSpEffBinlSplcngEffAry, @tmpJunctIntrvlCnsrvtnScrAry, @tmpSpEffBinCnsrvtnScrAry);

	my $totalJunctNum = keys %senseSplicingEfficiencySortingHsh;
	my $intervalJunctCounted = 0;
	my $totalJunctCounted = 0;
	my $intervalNum = 0;
	my %tmpConcateSeqHsh;
	my $allPdfPath = "";
	
	my %binEntropyHsh;
	
	foreach my $junctStr (sort {$senseSplicingEfficiencySortingHsh{$b} <=> $senseSplicingEfficiencySortingHsh{$a}} keys %senseSplicingEfficiencySortingHsh) {
		$intervalJunctCounted++;
		$totalJunctCounted++;
		
		print "Processing $totalJunctCounted of $totalJunctNum\r";
		
		my $leftFlankSeq = ${$NGSJunctInfoHsh{$junctStr}}{"leftFlankSeq"};
		my $rightFlankSeq = ${$NGSJunctInfoHsh{$junctStr}}{"rightFlankSeq"};
		my $oriLen = length $leftFlankSeq;
		my $trimLen = $oriLen-$trimExonEnd-$trimIntronEnd;
		my $leftExonSeq = substr $leftFlankSeq, 0, $boundWidth;
		my $rightExonSeq = substr $rightFlankSeq, $boundWidth+2;
		my $splicingEffType = ${$NGSJunctInfoHsh{$junctStr}}{"splicingEffType"};
		
		if ($splicingEffType eq "senseOften") {
			print EXONSEQFASTA ">$junctStr\n";
			print EXONSEQFASTA "$leftExonSeq"."$rightExonSeq\n";
		}
		
		die "the trimmed the length of the flank sequence is smaller than zero" if ($trimLen <= 0) ;
		my $trimLeftFlankSeq = substr $leftFlankSeq, $trimExonEnd, $trimLen;
		my $trimRightFlankSeq = substr $rightFlankSeq, $trimIntronEnd, $trimLen;
		my $concateSeq = $trimLeftFlankSeq.$trimRightFlankSeq;
		my $senseSplicingEfficiency = ${$NGSJunctInfoHsh{$junctStr}}{"senseSplicingEfficiency"};
		
		#---calculated the consensusValue;
		my $consensusValue = 0;
		my @concateSeqSplt = split //, $concateSeq;
		my $posSeq = "";
		my $minSum = 0;
		my $maxSum = 0;
		my $ntPrptnInRefPSSMSum = 0;
		
		foreach my $i (0..$#concateSeqSplt) {
			my $pos = $i + 1;
			if (exists $posFactorHsh{$pos}) {
				my $ntInConcateSeq = $concateSeqSplt[$i];
				$posSeq .= $ntInConcateSeq;
				my $ntPrptnInRefPSSM = ${$refPSSMHsh{$pos}}{$ntInConcateSeq};
				my $max = ${$refPSSMHsh{$pos}}{"max"};
				my $min = ${$refPSSMHsh{$pos}}{"min"};
				my $meanEntropy = ${$refPSSMHsh{$pos}}{"meanEntropy"};
				my $lowEntropy = ${$refPSSMHsh{$pos}}{"lowEntropy"};
				my $highEntropy = ${$refPSSMHsh{$pos}}{"highEntropy"};
				my $factor = ${$refPSSMHsh{$pos}}{"factor"};

				$minSum += $min * $factor;
				$maxSum += $max * $factor;
				$ntPrptnInRefPSSMSum += $ntPrptnInRefPSSM * $factor;

				#$minSum += $min;
				#$maxSum += $max;
				#$ntPrptnInRefPSSMSum += $ntPrptnInRefPSSM;

				#$minSum += $min * $meanEntropy;
				#$maxSum += $max * $meanEntropy;
				#$ntPrptnInRefPSSMSum += $ntPrptnInRefPSSM * $meanEntropy;

				#$minSum += $min * $meanEntropy * $meanEntropy;
				#$maxSum += $max * $meanEntropy * $meanEntropy;
				#$ntPrptnInRefPSSMSum += $ntPrptnInRefPSSM * $meanEntropy * $meanEntropy;

				#$consensusValue += (($ntPrptnInRefPSSM - $min)/($max - $min));
			}
		}
		$consensusValue = (($ntPrptnInRefPSSMSum - $minSum)/($maxSum - $minSum));
		
		${$NGSJunctInfoHsh{$junctStr}}{"consensusValue"} =  $consensusValue;
		
		#$consensusValue = ((($ntPrptnInRefPSSMSum*$ntPrptnInRefPSSMSum) - ($minSum*$minSum))/(($maxSum*$maxSum) - ($minSum*$minSum)));
		#$consensusValue = (($ntPrptnInRefPSSMSum - $minSum)/($maxSum - $minSum))*(($ntPrptnInRefPSSMSum - $minSum)/($maxSum - $minSum));
		#$consensusValue = (($ntPrptnInRefPSSMSum)/($maxSum));
		#$consensusValue = (($ntPrptnInRefPSSMSum)/($maxSum))*(($ntPrptnInRefPSSMSum)/($maxSum));

		my $logSenseSplicingEfficiency = log ($senseSplicingEfficiency)/log(10);
		print INDIVCNSCRLOG $logSenseSplicingEfficiency."\t".$consensusValue."\t".$senseSplicingEfficiency."\t".$posSeq."\n";
		
=pod
		my @similarityAry = ();
		foreach my $refPSSMJunctStr (keys %refPSSMPosSeqHsh) {
			my @refPosSeqSplt = split //, $refPSSMPosSeqHsh{$refPSSMJunctStr};
			my @posSeqSplt = split //, $posSeq;
			my $match = 0;
			for my $pos (0..$#refPosSeqSplt) {
				$match ++ if $refPosSeqSplt[$pos] eq $posSeqSplt[$pos];
			}
			my $similarity = 100*$match/@refPosSeqSplt;
			push @similarityAry, $similarity;
		}
		
		my ($meanSimilarity, $SD) = calculateStandardDeviationAndMean(\@similarityAry);
		print SMLRTYLOG $logSenseSplicingEfficiency."\t".$meanSimilarity."\t".$senseSplicingEfficiency."\t".$posSeq."\n";
=cut
		$tmpConcateSeqHsh{$junctStr} = $concateSeq;
		push @tmpJunctIntrvlSplcngEffAry, $senseSplicingEfficiency;
		push @tmpSpEffBinlSplcngEffAry, $senseSplicingEfficiency;
		push @tmpJunctIntrvlCnsrvtnScrAry, $consensusValue;
		push @tmpSpEffBinCnsrvtnScrAry, $consensusValue;

		if (($logSenseSplicingEfficiency <= $lowerbinLimit) or ($totalJunctCounted == $totalJunctNum)) {
			my $binMid = $lowerbinLimit - ($binWidthLogSenseSplicingEfficiency/2);
			my ($meanSplicigEff, $SDSplicigEff) = calculateStandardDeviationAndMean(\@tmpSpEffBinlSplcngEffAry);
			my ($meanCnsvtnScr, $SDCnsvtnScr) = calculateStandardDeviationAndMean(\@tmpSpEffBinCnsrvtnScrAry);
			my $logMeanSplicigEff = log($meanSplicigEff)/log(10);
			my $indivCScoreStr = join "\t", @tmpSpEffBinCnsrvtnScrAry;
			my $n = @tmpSpEffBinCnsrvtnScrAry;
			my $binTag = $upperbinLimit."To".$lowerbinLimit."_n=$n";
			print EFFBINCNSCRLOG $binTag."\t".$logMeanSplicigEff."\t".$meanCnsvtnScr."\t".$SDCnsvtnScr."\t".$n."\n";

			my $binSeqFastaFilePath = "$outDir/splicingEffInterval/fasta/$binTag.fasta";
			open (BINSEQ, ">$binSeqFastaFilePath");
			foreach my $binJunctStr (keys %tmpConcateSeqHsh) {
				my $binConcateSeq = $tmpConcateSeqHsh{$binJunctStr};
				print BINSEQ ">$binJunctStr\n";
				print BINSEQ "$binConcateSeq\n";
			}
			close BINSEQ;

			my $entropyPath = "$outDir/splicingEffInterval/entropy/$binTag.entropy.txt";
			my $seqPath = $binSeqFastaFilePath;
			my $logoPath = "$outDir/splicingEffInterval/pdf/$binTag.pdf";
		
			$allPdfPath .= " \'$logoPath\'";
		
			my $title = $binTag;
			my $numSeq = `sed -n \'\$=\' $seqPath`;
			$numSeq = $numSeq/2;
			$title .= "_n=$numSeq";
			my $weblogoEntropyCMD = "weblogo -f $seqPath -F logodata -n $length -s large --title $title -A rna -c classic >$entropyPath";
			my $grepCMD = "ps -ef | grep weblogo | grep -v grep | grep -v perl";
			runAndCheckSerialTask($grepCMD, "weblogo", $weblogoEntropyCMD, "$outDir/error.log.txt");

			$weblogoEntropyCMD = "weblogo -f $seqPath -n $length -F pdf -s large --title $title -A rna -c classic >$logoPath";
			$grepCMD = "ps -ef | grep weblogo | grep -v grep | grep -v perl";
			runAndCheckSerialTask($grepCMD, "weblogo", $weblogoEntropyCMD, "$outDir/error.log.txt");
			
			open (INTENTROPY, "$entropyPath");
			while (my $theLine = <INTENTROPY>) {
				chomp $theLine;
				$theLine =~ s/ //g;
				next if (($theLine =~ m/^#/) or (length $theLine < 2));
				my ($pos, $ACount, $CCount, $GCount, $UCount, $meanEntropy, $lowEntropy, $highEntropy, $weight) = split /\t/, $theLine;
				${$binEntropyHsh{$pos}}{$logMeanSplicigEff} = $meanEntropy;
			}
			close INTENTROPY;

			$upperbinLimit = $lowerbinLimit;
			$lowerbinLimit = $lowerbinLimit - $binWidthLogSenseSplicingEfficiency;
			@tmpSpEffBinCnsrvtnScrAry = ();
			@tmpSpEffBinlSplcngEffAry = ();
			%tmpConcateSeqHsh = ();
		}
		
		if ($intervalJunctCounted >= $numJunctInterval) {
			$intervalNum++;
			my ($meanSplicigEff, $SDSplicigEff) = calculateStandardDeviationAndMean(\@tmpJunctIntrvlSplcngEffAry);
			my ($meanCnsvtnScr, $SDCnsvtnScr) = calculateStandardDeviationAndMean(\@tmpJunctIntrvlCnsrvtnScrAry);
			my $binTag = "int_$intervalNum.mean$meanSplicigEff";

			my $logMeanSplicigEff = log ($meanSplicigEff)/log(10);
			my $indivCScoreStr = join "\t", @tmpJunctIntrvlCnsrvtnScrAry;
			print INTRVCNSCRLOG $binTag."\t".$logMeanSplicigEff."\t".$indivCScoreStr."\n";
			
			#---reset the counter
			$intervalJunctCounted = 0;
			@tmpJunctIntrvlSplcngEffAry = ();
			@tmpJunctIntrvlCnsrvtnScrAry = ();
		}
	}
	
	close EXONSEQFASTA;
	close INDIVCNSCRLOG;
	close INTRVCNSCRLOG;
	close SMLRTYLOG;
	close EFFBINCNSCRLOG;

	my $allBinEntropyPath = "$outDir/splicingEffInterval/log/allBinEntropy.log.txt";
	open (ALLBINENTHROPY, ">$allBinEntropyPath");
	foreach my $pos (sort {$a <=> $b} keys %binEntropyHsh) {
		print ALLBINENTHROPY "pos";
		foreach my $logMeanSplicigEff (sort {$a <=> $b} keys %{$binEntropyHsh{$pos}}) {
			print ALLBINENTHROPY "\t".$logMeanSplicigEff;
		}
		print ALLBINENTHROPY "\n";
		last;
	}
	foreach my $pos (sort {$a <=> $b} keys %binEntropyHsh) {
		print ALLBINENTHROPY $pos;
		foreach my $logMeanSplicigEff (sort {$a <=> $b} keys %{$binEntropyHsh{$pos}}) {
			my $meanEntropy = ${$binEntropyHsh{$pos}}{$logMeanSplicigEff};
			print ALLBINENTHROPY "\t".$meanEntropy;
		}
		print ALLBINENTHROPY "\n";
	}
	close ALLBINENTHROPY;
	
	print "\nMerging all pdfs.\n";
	my $mergePdfCMD = "gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=\'$outDir/splicingEffInterval/pdf/all.merged.pdf\' $allPdfPath";
	my $grepCMD = "ps -ef | grep $allPdfPath | grep -v grep";
	runAndCheckSerialTask($grepCMD, $mergePdfCMD, $mergePdfCMD, "$outDir/error.log.txt");
	
	my $dirForDremePath = "$outDir/dreme/shuffle/";
	system "mkdir -p -m 777 $dirForDremePath";
	my $dremeCMD = "dreme -oc $dirForDremePath -p $exonSeqFastaPath -maxk 10 -mink 4 -e 0.001";
	$grepCMD = "ps -ef | grep dreme | grep -v grep | grep -v perl";
	my $grepStr = "dreme -oc $dirForDremePath";
	runAndCheckSerialTask($grepCMD, $grepStr, $dremeCMD, "$outDir/error.log.txt");
	plotDremeMotif($dirForDremePath, $dirForDremePath."./dreme.xml", $exonSeqFastaPath, "exonSeqOftenSplice", "exonSeqOftenSplice fullShuffleSearch");

	$dirForDremePath = "$outDir/dreme/neg/";
	system "mkdir -p -m 777 $dirForDremePath";
	my $randCDSSeqPath = "$outDir/dreme/fasta/randCDS.fasta";
	my $randSeqLen = 100;
	my $randSeqNum = 2000;
	randomGenDNASeq(\%intronlessGeneCodingSeqHsh, $randSeqLen, $randSeqNum, $randCDSSeqPath);

	$dremeCMD = "dreme -oc $dirForDremePath -p $exonSeqFastaPath -n $randCDSSeqPath -maxk 10 -mink 4 -e 0.001";
	$grepCMD = "ps -ef | grep dreme | grep -v grep | grep -v perl";
	$grepStr = "dreme -oc $dirForDremePath";
	runAndCheckSerialTask($grepCMD, $grepStr, $dremeCMD, "$outDir/error.log.txt");
	plotDremeMotif($dirForDremePath, $dirForDremePath."./dreme.xml", $exonSeqFastaPath, "exonSeqOftenSplice", "exonSeqOftenSplice vs random intronless CDS");
	
	return \%NGSJunctInfoHsh;

}
########################################################################## plotDremeMotif
sub plotDremeMotif {
	
	#---subroutine dependency: getDremeMotif, readMultiFasta, GNUPlotMultipleColumnXYLines, getSeqComposition
	#---in/out: plotDremeMotif($dremeXMLToRead, $fastaSeq);
	
	my $subOutDir = $_[0];
	my $dremeXMLToRead = $_[1];
	my $fastaSeq = $_[2];
	my $fileTag = $_[3];
	my $title = $_[4];
	
	system "mkdir -p -m 777 $subOutDir/";

	my $motifInfoHsh_ref = getDremeMotif($dremeXMLToRead);
	my %motifInfoHsh = %{$motifInfoHsh_ref};
	if ((keys %motifInfoHsh) >= 1) {
		my $fastaHsh_ref = readMultiFasta($fastaSeq);
		my $motifOccurenceFilePath = $subOutDir."/$fileTag.motifOccurence.txt";
		my $motifOccurencePlotFilePath = $subOutDir."/$fileTag.motifOccurence.pdf";
		scanForMotifOccurence($motifInfoHsh_ref, $fastaHsh_ref, $motifOccurenceFilePath);
		GNUPlotMultipleColumnXYLines($motifOccurencePlotFilePath, $motifOccurenceFilePath, $title, "percentage");
	} else {
		open TMP, ">$subOutDir/$fileTag.has.no.significant.motifs.txt"; close TMP;
	}
	
}
########################################################################## randomGenDNASeq
sub getDremeMotif {

	my $dremeXMLToRead = $_[0];
	
	my %motifInfoHsh;
	my ($motifSeq, $evalue);
	open (DREMEXML, "$dremeXMLToRead");
	while (my $theLine = <DREMEXML>) {
		chomp $theLine;
		if ($theLine =~ m/<motif id=/) {
			my @theLineSplt = split / |\>/, $theLine;
			foreach my $arg (@theLineSplt) {
				$arg =~ s/\"//g;
				if ($arg =~ m/^seq=/) {$motifSeq = substr ($arg, index ($arg, "=")+1);}
				elsif ($arg =~ m/^evalue=/) {$evalue = substr ($arg, index ($arg, "=")+1);}
			}
			${$motifInfoHsh{$motifSeq}}{"evalue"} = $evalue;
		}
		if ($theLine =~ m/<match seq=/) {
			my @theLineSplt = split / |\>/, $theLine;
			foreach my $arg (@theLineSplt) {
				if ($arg =~ m/^seq=/) {
					$arg =~ s/\"//g;
					my $matchSeq = substr ($arg, index ($arg, "=")+1); 
					push @{${$motifInfoHsh{$motifSeq}}{"matchSeq"}}, $matchSeq;
				}
			}
		}
	}
	close (DREMEXML);
	
	my $motifNum = keys %motifInfoHsh;
	
	print "$motifNum motif stored.\n";
	
	return \%motifInfoHsh;

}
########################################################################## scanForMotifOccurence
sub scanForMotifOccurence {
	
	#---subroutine dependency: 
	#---in/out: scanForMotifOccurence($var);
	
	my %motifInfoHsh = %{$_[0]}; #---push @{${$motifInfoHsh{$motifSeq}}{"matchSeq"}}, $matchSeq;
	my %fastaHsh = %{$_[1]};
	my $outFilePath = $_[2];
	
	my %motifHitHsh;
	my %seqHitInfoHsh;
	my $totalSeqNum = keys %fastaHsh;
	my %allPosHsh; #---just to store the position;

	foreach my $motifSeq (sort {${$motifInfoHsh{$a}}{"evalue"} <=> ${$motifInfoHsh{$b}}{"evalue"}} keys %motifInfoHsh) {
		print "Scanning motif occurence for $motifSeq\r";
		%allPosHsh = ();
		foreach my $seqName (keys %fastaHsh) {
			my %seqHitPosHsh = ();
			my $seq = $fastaHsh{$seqName};
			$seq =~ tr/Uu/Tt/;
			
			for my $pos (0..((length $seq)-1)) {
				${$seqHitPosHsh{$seqName}}{$pos} = 0;
				$allPosHsh{$pos}++;
				${${$motifHitHsh{$motifSeq}}{$pos}}{"count"} = 0 if (not exists ${${$motifHitHsh{$motifSeq}}{$pos}}{"count"}); #---build the motifHitHsh base
			}

			foreach my $matchSeq (@{${$motifInfoHsh{$motifSeq}}{"matchSeq"}}) {
				my $motifLen = length $motifSeq;
				my @startPosAry;
				my $offset = 0;
				my $startPos = index(uc ($seq), $matchSeq, $offset);#---UC for upper case as dreme always output uppercase
				while ($startPos != -1) {
					push @startPosAry, $startPos;
					$offset = $startPos + 1;
					$startPos = index(uc ($seq), $matchSeq, $offset);
				}
				
				foreach my $startPos (@startPosAry) {
					foreach my $hitPos ($startPos..($startPos+$motifLen-1)) {
						${$seqHitPosHsh{$seqName}}{$hitPos}++;
					}
				}
			} #---end foreach my $matchSeq (@{${$motifInfoHsh{$motifSeq}}{"matchSeq"}}) {
			
			foreach my $hitPos (keys %{$seqHitPosHsh{$seqName}}) {
				${${$motifHitHsh{$motifSeq}}{$hitPos}}{"count"}++ if (${$seqHitPosHsh{$seqName}}{$hitPos} > 0);
			}
		}#--- end of foreach my $seqName (keys %fastaHsh) {
		
		foreach my $pos (keys %{$motifHitHsh{$motifSeq}}) {
			my $count = ${${$motifHitHsh{$motifSeq}}{$pos}}{"count"};
			my $pct = sprintf "%.06f", 100*$count/$allPosHsh{$pos};
			${${$motifHitHsh{$motifSeq}}{$pos}}{"pct"} = $pct;
		}
	}#---end of foreach my $motifSeq (sort {${$motifInfoHsh{$a}}{"evalue"} <=> ${$motifInfoHsh{$b}}{"evalue"}} keys %motifInfoHsh) {
	print "......finished\n";

	open (MOTIFOCC, ">$outFilePath");
	print MOTIFOCC "pos";
	foreach my $motifSeq (sort {${$motifInfoHsh{$a}}{"evalue"} <=> ${$motifInfoHsh{$b}}{"evalue"}} keys %motifInfoHsh) {
		print MOTIFOCC "\t".$motifSeq."_".${$motifInfoHsh{$motifSeq}}{"evalue"};
	}
	print MOTIFOCC "\n";
	foreach my $pos (sort {$a <=> $b} keys %allPosHsh) {
		print MOTIFOCC $pos;
		foreach my $motifSeq (sort {${$motifInfoHsh{$a}}{"evalue"} <=> ${$motifInfoHsh{$b}}{"evalue"}} keys %motifInfoHsh) {
			print MOTIFOCC "\t".${${$motifHitHsh{$motifSeq}}{$pos}}{"pct"};
		}
		print MOTIFOCC "\n";
	}
	close MOTIFOCC;

	return;
}
########################################################################## GNUPlotMultipleColumnXYLines
sub GNUPlotMultipleColumnXYLines {

	#---GNUPlotMultipleColumnXYLines($plotFilePath, $plotDataPath, $title, $ylabel);
	#The file must look like this
	#X	Y1	Y2	Y3
	#1	45	872	68
	#2	45	87	53
	#3	68	97	1
	#4	10	2	60
	
	my $plotFilePath = $_[0];
	my $plotDataPath = $_[1];
	my $title = $_[2];
	my $ylabel = $_[3];
	my $extraCmd = $_[4];
	
	$extraCmd = "" if (not defined $extraCmd);
	
	$plotFilePath .= ".pdf" if ($plotFilePath !~ m/\.pdf$/);

	my @filePathSplt = split /\//, $plotFilePath;
	my $fileName = $filePathSplt[-1];

	my @allPlotCmdAry;
	my @indivPlotCmdAry;
	my $xlabel;
	open PLOTDATA, "$plotDataPath" ;
	while (my $theLine = <PLOTDATA>) {
		chomp $theLine;
		my @theLineSplt = split /\t/, $theLine;
		$xlabel = $theLineSplt[0];
		foreach my $i (1..$#theLineSplt) {
			my $legend = $theLineSplt[$i];
			my $YColumn = $i+1;
			push @allPlotCmdAry, "\'$plotDataPath\' using 1:$YColumn with lines title \"$legend\" ";
			push @indivPlotCmdAry, "plot \'$plotDataPath\' using 1:$YColumn with lines title \"$legend\"; ";
		}
		last;
	}
	close PLOTDATA;
	
	my $allPlotCmdstr = join ",", @allPlotCmdAry;
	my $indivPlotCmdStr = join " ", @indivPlotCmdAry;

	$allPlotCmdstr = "plot ".$allPlotCmdstr.";";
	
	print "Running GNUPlotMultipleColumnXYLines for $fileName.\n";
	
	#---do the GNUPLOT
	open (GNUPLOT, "|gnuplot");
	print GNUPLOT <<EOPLOT;
	set terminal postscript color solid
	set output "| ps2pdf - $plotFilePath 2>/dev/null";
	unset logscale x;
	unset logscale y;
	$extraCmd;
	set xlabel "$xlabel";
	set ylabel "$ylabel";
	set title "$title";
	$allPlotCmdstr;
	$indivPlotCmdStr;
EOPLOT
	close(GNUPLOT);
}
########################################################################## genRandomCodingSequence
sub getCodingSequence {
	
	#---my $intronlessGeneCodingSeqHsh_ref = getCodingSequence(\%cntgSeqHsh, \%strndByGeneHsh, \%exonRngByGeneHsh, \%cntgByGeneHsh, \%geneCtgryHsh);
	
	my %cntgSeqHsh = %{$_[0]}; 
	my %strndByGeneHsh = %{$_[1]}; #--- $strndByGeneHsh{$geneID} = $geneStrd;
	my %exonRngByGeneHsh = %{$_[2]}; #---  ${${$exonRngByGeneHsh{$geneID}}{$exonID}}{"start"} = $featureStart;
	my %cntgByGeneHsh = %{$_[3]}; #--- $cntgByGeneHsh{$geneID} = $seq;
	my %geneCtgryHsh = %{$_[4]};

	my %intronlessGeneCodingSeqHsh;
	foreach my $geneID (keys %exonRngByGeneHsh) {
	
		my $cntg = $cntgByGeneHsh{$geneID};
		my $cntgSeq = $cntgSeqHsh{$cntg};
		my $geneStrd = $strndByGeneHsh{$geneID};
			
		#---get the boundaries
		my (@tmpExonBoundAry);
		my @originalExonBoundAry;
		my $tmpExonStr = "";
		my $allExonStart = "";
		my $allExonEnd = "";
		my $exonNum = keys %{$exonRngByGeneHsh{$geneID}};
		
		if ($exonNum == 1){#---no intron
			my $seqToReturn = "";
			foreach my $exonID (keys %{$exonRngByGeneHsh{$geneID}}) {
				push  @tmpExonBoundAry, ${${$exonRngByGeneHsh{$geneID}}{$exonID}}{"start"};
				push  @tmpExonBoundAry, ${${$exonRngByGeneHsh{$geneID}}{$exonID}}{"end"};
			}
			for (my $i = 0; $i < $#tmpExonBoundAry; $i += 2) {#---$i = 0, 2, 4, 6 is there're are 8 elements
				my $exonStart = $tmpExonBoundAry[$i];
				my $exonEnd = $tmpExonBoundAry[$i+1];
				my $length = $exonEnd-$exonStart+1;
				my $offset = $exonStart-1;
				$offset = 0 if ($offset < 0);
				$seqToReturn .= substr ($cntgSeq, $offset, $length);
			}
			
			if ($geneStrd eq "-") {
				$seqToReturn = reverse $seqToReturn;
				$seqToReturn =~ tr/ACGTacgt/TGCAtgca/;
			}
			
			$intronlessGeneCodingSeqHsh{$geneID} = $seqToReturn;
			#print $geneID."\n".$seqToReturn."\n";
		}
	}#--- end of foreach my $geneID (sort {$a cmp $b} keys %exonRngByGeneHsh) {
	
	return \%intronlessGeneCodingSeqHsh;
	
}
########################################################################## randomGenDNASeq
sub randomGenDNASeq {
	
	#---my $randSeqHsh_ref = randomGenDNASeq(\%refSeqHsh, $seqLen, $seqNum, $fastaPath);
	my %refSeqHsh = %{$_[0]};
	my $seqLen = $_[1];
	my $seqNum = $_[2];
	my $fastaPath = $_[3];#---input the path for printing or use "no" to swicth off;
	
	print "Generating $seqNum random sequences of $seqLen in length.\n";
	
	open (FASTA, ">$fastaPath") if ($fastaPath ne "no");
	
	my %randSeqHsh;
	my $conCatSeq = "";
	
	foreach my $seqName (keys %refSeqHsh) {
		my $seq = $refSeqHsh{$seqName};
		$conCatSeq .= $seq;
	}
	
	my $totalConCatLength = length $conCatSeq;
	
	for my $round (1..$seqNum) {
		my $validSeq = "no";

		while ($validSeq eq "no") {#---loop until not out of ctng rng

			my $seqStartPos = int (rand $totalConCatLength); #---rand from 0 to $cntgLen-1;
			next if (($seqStartPos + $seqLen) > $totalConCatLength);
			my $extractSeq = substr $conCatSeq, $seqStartPos, $seqLen;
		 	next if ($extractSeq =~ m/[^ATGCatgc]/);
		 	$validSeq = "yes";
		 	
		 	$randSeqHsh{"randSeq".$round} = $extractSeq;
		 	
		 	if ($fastaPath ne "no") {
			 	print FASTA ">"."randSeq".$round."\n";
			 	print FASTA $extractSeq."\n";
			}
		}
	}
	
	return \%randSeqHsh;
}
########################################################################## printGeneBasedJunctCount
sub printGeneBasedJunctCount {

	#----printGeneBasedJunctCount($geneBasedJunctCountHsh_ref, $refNameByGeneHsh_ref);
	
	my %geneBasedJunctCountHsh = %{$_[0]};
	my %refNameByGeneHsh = %{$_[1]};
	
	print "Printing gene based junction count\n";
	open (GENEJCOUNT, ">$outDir/geneBasedInfo/geneBasedJunctCount.log.txt");
	print GENEJCOUNT join "", ((join "\t", ("geneID", "geneName","senseTotal", "senseOften", "senseRare", "antisenseTotal", "antisenseOften", "antisenseRare")), "\n");
	foreach my $geneID (sort {$a cmp $b} keys %refNameByGeneHsh) {
		my $geneName = $refNameByGeneHsh{$geneID};
		print GENEJCOUNT $geneID."\t".$geneName;
		foreach my $junctType (("senseTotal", "senseOften", "senseRare", "antisenseTotal", "antisenseOften", "antisenseRare")) {
			my $count = 0;
			$count = ${$geneBasedJunctCountHsh{$geneID}}{$junctType} if exists ${$geneBasedJunctCountHsh{$geneID}}{$junctType};
			print GENEJCOUNT "\t".$count;
		}
		print GENEJCOUNT "\n";
	}
	close GENEJCOUNT;
}
########################################################################## calculateProportionOfCodingAtSplicingEfficienyInterval
sub calculateProportionOfCodingAtSplicingEfficienyInterval {
	
	#---calculateProportionOfCodingAtSplicingEfficienyInterval($NGSJunctInfoHsh_ref, $numJunctInterval, $minReadNumForPrprtnORF);
	
	my %NGSJunctInfoHsh = %{$_[0]};
	my $numJunctInterval = $_[1];
	my $minReadNumForPrprtnORF = $_[2];
	
	my (@splcingEffAry, %altIsofmFullORFHsh);
	my %altJunctSpEffHsh;#---for sorting purpose since senseSplicingEfficiency can sometimes be non integer;
	foreach my $junctStr (keys %NGSJunctInfoHsh) {
		my $refOrAlt = ${$NGSJunctInfoHsh{$junctStr}}{"refOrAlt"};
		my $readNum = ${$NGSJunctInfoHsh{$junctStr}}{"readNum"};
		next if (($readNum < $minReadNumForPrprtnORF) or ($refOrAlt ne "alt"));
		my $senseSplicingEfficiency = ${$NGSJunctInfoHsh{$junctStr}}{"senseSplicingEfficiency"};
		$altJunctSpEffHsh{$junctStr} = $senseSplicingEfficiency;
	}

	open (SPEFFORF, ">$outDir/prprtnORFSpEffIntrvl.log.txt");
	my $numJunctCounted = 0;
	foreach my $junctStr (sort {$altJunctSpEffHsh{$b} <=> $altJunctSpEffHsh{$a}} keys %altJunctSpEffHsh) {
		$numJunctCounted++;
		my $senseSplicingEfficiency = ${$NGSJunctInfoHsh{$junctStr}}{"senseSplicingEfficiency"};
		my $altIsofmFullORF = ${$NGSJunctInfoHsh{$junctStr}}{"altIsofmFullORF"};
		my $altSplcType = ${$NGSJunctInfoHsh{$junctStr}}{"altSplcType"};
		push @splcingEffAry, $senseSplicingEfficiency;
		$altIsofmFullORFHsh{$altIsofmFullORF}++;
		
		if ($numJunctCounted == $numJunctInterval) {
			my ($meanSplcingEff, $SDSplcingEff) = calculateStandardDeviationAndMean(\@splcingEffAry);
			my $fullORFPrptn = $altIsofmFullORFHsh{"y"}/@splcingEffAry;
			print SPEFFORF $fullORFPrptn."\t".$meanSplcingEff."\t".$SDSplcingEff."\t".@splcingEffAry."\n";
			
			@splcingEffAry = ();
			%altIsofmFullORFHsh = ();
			$numJunctCounted = 0;
		}
	}
	close SPEFFORF;
}
########################################################################## calculateMeanJunctionNumWithinClusterRepJunctRdNumInterval
sub calculateMeanJunctionNumWithinClusterRepJunctRdNumInterval {

	#---calculateMeanJunctionNumWithinClusterRepJunctRdNumInterval($jClusterInfoHsh_ref, $intervalBin, $logScale);
	
	my %jClusterInfoHsh = %{$_[0]};
	my $intervalBin = $_[1];
	my $logScale = $_[2]; #----10 or 2, or any integer
	my %NGSJunctInfoHsh = %{$_[3]};
	
	my $intervalCutoff = $intervalBin;
	my $clusterProc = 0;
	my $totalClusterNum = keys %jClusterInfoHsh;
	my (@tmpRepRdNumAry, @tmpJunctNumAry, @tmpAltRdNumAry);
	
	my $likelyAllNoiseMaxCV = 0.8;
	
	open (CLRDNUMJNUM, ">$outDir/correlation/withinClusterRdNumIntvlVsJunctNum.log.txt");
	open (CLREPRDALTRD, ">$outDir/correlation/withinClusterRepRdNumIntvlVsAltRdNum.log.txt");
	open (CLRACTUALEPRDALTRD, ">$outDir/correlation/withinClusterWithAltJunctRepAltJunctRdRatioCutoff$likelyAllNoiseMaxCV.log.txt");
	open (INDIVACTUALEPRDALTRD, ">$outDir/correlation/individualAltJunctRepAltJunctRdRatio.log.txt");

	my %repAltJunctRdRatioHsh;
	my @repAltRdNumInClusterRatioAry;

	foreach my $clusterName (sort {${$jClusterInfoHsh{$a}}{"mostSupportReadNum"} <=> ${$jClusterInfoHsh{$b}}{"mostSupportReadNum"}} keys %jClusterInfoHsh) {
		$clusterProc++;

		my $ovrlpCtrgy = ${$jClusterInfoHsh{$clusterName}}{"ovrlpCtrgy"};
		my $ovrlpSense = ${$jClusterInfoHsh{$clusterName}}{"ovrlpSense"};
		my $ovrlpAntisense = ${$jClusterInfoHsh{$clusterName}}{"ovrlpAntisense"};
		my $superCtgry = ${$jClusterInfoHsh{$clusterName}}{"superCtgry"};
		my $splicingRatio = ${$jClusterInfoHsh{$clusterName}}{"splicingRatio"};
		my $splicingEffType = ${$jClusterInfoHsh{$clusterName}}{"splicingEffType"};
		my $hitRef = ${$jClusterInfoHsh{$clusterName}}{"hitRef"};
		my $unq = ${$jClusterInfoHsh{$clusterName}}{"unq"};
		my $majorJunctNum = ${$jClusterInfoHsh{$clusterName}}{"majorJunctNum"};
		my $totalJunctNum = ${$jClusterInfoHsh{$clusterName}}{"totalJunctNum"};
		my $majorCluster = ${$jClusterInfoHsh{$clusterName}}{"majorCluster"};
		my $prominentJunct = ${$jClusterInfoHsh{$clusterName}}{"prominentJunct"};
		my $mostSupportReadNum = ${$jClusterInfoHsh{$clusterName}}{"mostSupportReadNum"};
		my $superCluster = ${$jClusterInfoHsh{$clusterName}}{"superCluster"};
		my $totalRdNumInCluster = ${$jClusterInfoHsh{$clusterName}}{"totalRdNumInCluster"};
		my $altRdNumInCluster = $totalRdNumInCluster - $mostSupportReadNum;
		
		if ($ovrlpSense ne "null") {#--- on sense mRNA only
			if ($majorJunctNum >= 1) {#----constitutive clusters
				my $logMostSupportReadNum = log($mostSupportReadNum)/log($logScale);
				my $likelyAllNoise = "yes";
				foreach my $altJunctStr (@{${$jClusterInfoHsh{$clusterName}}{"allJunctStr"}}) {
					next if ($altJunctStr eq $prominentJunct);
					my $altJunctRdNum = ${$NGSJunctInfoHsh{$altJunctStr}}{"readNum"};
					my $logAltJunctRdNum = log($altJunctRdNum)/log($logScale);
					my $logRepAltJunctRdRatio = $logAltJunctRdNum/$logMostSupportReadNum;
					${$repAltJunctRdRatioHsh{$altJunctStr}} = $logRepAltJunctRdRatio;
					my $altJunctRdPct = $altJunctRdNum/$totalRdNumInCluster;
					my $repAltJunctRdRatio = $altJunctRdNum/$mostSupportReadNum;
					print INDIVACTUALEPRDALTRD $altJunctStr."\t".$altJunctRdNum."\t".$mostSupportReadNum."\t".$totalRdNumInCluster."\t".$logRepAltJunctRdRatio."\t".$repAltJunctRdRatio."\t".$altJunctRdPct."\n";

					my $consensusValue = ${$NGSJunctInfoHsh{$altJunctStr}}{"consensusValue"};
					if (($consensusValue > $likelyAllNoiseMaxCV) or ($consensusValue eq "null")) {#---$consensusValue eq "null" might equalt to null if the alt junct is not unique depends on prior definitions
						$likelyAllNoise = "no";
					}
				}
				
				my $logAltRdNumInCluster = 0;
				$logAltRdNumInCluster = log($altRdNumInCluster)/log($logScale) if ($altRdNumInCluster > 0);
				push @tmpRepRdNumAry, $logMostSupportReadNum;
				push @tmpJunctNumAry, $totalJunctNum;
				push @tmpAltRdNumAry, $logAltRdNumInCluster;
				
				if ($logAltRdNumInCluster > 0) {
					my $logRepAltRdNumInClusterRatio = $logAltRdNumInCluster/$logMostSupportReadNum;
					my $repAltRdNumInClusterRatio = $altRdNumInCluster/$mostSupportReadNum;
					push @repAltRdNumInClusterRatioAry, $logRepAltRdNumInClusterRatio;
					print CLRACTUALEPRDALTRD $likelyAllNoise."\t".$mostSupportReadNum."\t".$altRdNumInCluster."\t".$logRepAltRdNumInClusterRatio."\t".$repAltRdNumInClusterRatio."\n";
				}
				
				if (($logMostSupportReadNum >= $intervalCutoff) or ($clusterProc == $totalClusterNum)) {
					my ($meanRepRdNum, $SDRepRdNum) = calculateStandardDeviationAndMean(\@tmpRepRdNumAry);
					my ($meanJunctNum, $SDJunctNum) = calculateStandardDeviationAndMean(\@tmpJunctNumAry);
					my ($meanAltRdNum, $SDAltRdNum) = calculateStandardDeviationAndMean(\@tmpAltRdNumAry);
					print CLRDNUMJNUM $meanRepRdNum."\t".$meanJunctNum."\t".$SDJunctNum."\t".@tmpJunctNumAry."\n";
					print CLREPRDALTRD $meanRepRdNum."\t".$meanAltRdNum."\t".$SDAltRdNum."\t".@tmpAltRdNumAry."\n";
					@tmpRepRdNumAry = ();
					@tmpJunctNumAry = ();
					@tmpAltRdNumAry = ();
					$intervalCutoff += $intervalBin;
				}
			}
		}
	}
	
	close CLRDNUMJNUM;
	close CLREPRDALTRD;
	close CLRACTUALEPRDALTRD;
	close INDIVACTUALEPRDALTRD;
}
########################################################################## calculateProportionOfIntronCreationAtAbundanceInterval
sub calculateProportionOfIntronCreationAtAbundanceInterval  {

	#---calculateProportionOfIntronCreationAtAbundanceInterval($NGSJunctInfoHsh_ref, $countFturCovInfoHsh_ref, $refStrndHsh_ref, $intervalBin, $logScale);

	my %NGSJunctInfoHsh = %{$_[0]};
	my %countFturCovInfoHsh = %{$_[1]};
	my %refStrndHsh = %{$_[2]};
	my $intervalBin = $_[3];
	my $logScale = $_[4]; #----10 or 2, or any integer
	
	my %geneIntronCreationInfoHsh;

	print "Checking intron creation and transcript abundance correlation";
	open (INTNCRTABNRDNUM, ">$outDir/correlation/intronCreationRDNumAtAbundanceInterval.log.txt");
	open (INTNCRTABNJNUM, ">$outDir/correlation/intronCreationJunctNumAtAbundanceInterval.log.txt");

	foreach my $junctStr (keys %NGSJunctInfoHsh) {
		my $readNum = ${$NGSJunctInfoHsh{$junctStr}}{"readNum"};
		my $altSplcType = ${$NGSJunctInfoHsh{$junctStr}}{"altSplcType"};
		my $onRefGene = ${$NGSJunctInfoHsh{$junctStr}}{"onRefGene"};
		
		if ($altSplcType eq "intronCreate") {
			${$geneIntronCreationInfoHsh{$onRefGene}}{"intronCreateRdNum"} = 0 if not exists $geneIntronCreationInfoHsh{$onRefGene};
			${$geneIntronCreationInfoHsh{$onRefGene}}{"intronCreateJNum"}++;
			${$geneIntronCreationInfoHsh{$onRefGene}}{"intronCreateRdNum"} += $readNum;
		}
	}
	
	my $intervalCutoff = $intervalBin;
	my $procGene = 0;
	my $totalGeneNum = keys %refStrndHsh;
	my (@tmpICRdNumAry, @tmpICJunctNumAry, @tmpGeneCovNtAry);
	
	foreach my $geneID (sort {${$countFturCovInfoHsh{$a}}{"bothCov"} <=> ${$countFturCovInfoHsh{$b}}{"bothCov"}} keys %countFturCovInfoHsh) {
		$totalGeneNum++;
		my $geneCovNt = ${$countFturCovInfoHsh{$geneID}}{"bothCov"}/${$countFturCovInfoHsh{$geneID}}{"length"};
		next if ((not exists $refStrndHsh{$geneID}) or ($geneCovNt == 0));
		my $logGeneCovNt = log($geneCovNt)/log($logScale);
		if ($totalGeneNum == 1) {
			while ($intervalCutoff < $logGeneCovNt) {
				$intervalCutoff += $intervalBin;
			}
		}
		my $logIntronCreateRdNum = 0;
		my $intronCreateJNum = 0;
		if (exists $geneIntronCreationInfoHsh{$geneID}) {
			$intronCreateJNum = ${$geneIntronCreationInfoHsh{$geneID}}{"intronCreateJNum"};
			my $intronCreateRdNum = ${$geneIntronCreationInfoHsh{$geneID}}{"intronCreateRdNum"};
			$logIntronCreateRdNum = log($intronCreateRdNum)/log($logScale);
		}
		
		push @tmpICRdNumAry, $logIntronCreateRdNum;
		push @tmpICJunctNumAry, $intronCreateJNum;
		push @tmpGeneCovNtAry, $logGeneCovNt;
		
		if (($logGeneCovNt > $intervalCutoff) or ($procGene == $totalGeneNum)) {
			my ($meanICRdNum, $SDICRdNum) = calculateStandardDeviationAndMean(\@tmpICRdNumAry);
			my ($meanICJunctNum, $SDICJunctNum) = calculateStandardDeviationAndMean(\@tmpICJunctNumAry);
			my ($meanGeneCovNt, $SDGeneCovNt) = calculateStandardDeviationAndMean(\@tmpGeneCovNtAry);
			print INTNCRTABNRDNUM $meanGeneCovNt."\t".$meanICRdNum."\t".$SDICRdNum."\t".@tmpICRdNumAry."\n";
			print INTNCRTABNJNUM $meanGeneCovNt."\t".$meanICJunctNum."\t".$SDICJunctNum."\t".@tmpICJunctNumAry."\n";
			@tmpICRdNumAry = ();
			@tmpICJunctNumAry = ();
			@tmpGeneCovNtAry = ();
			$intervalCutoff += $intervalBin;
		}
	}
	
	close INTNCRTABNRDNUM;
	close INTNCRTABNJNUM;
}
########################################################################## plotSplicingEfficiencyVersusSiteEntropy
sub dremeNonStochasticJunctions {
	
	#--- dremeNonStochasticJunctions(\%NGSJunctInfoHsh, $intronlessGeneCodingSeqHsh_ref, $minPrmntReadNumRatioAsNonStochasticAltSite, $minSenseSplicingEfficiencyAsNonStochasticIntronCreate);
	
	my %NGSJunctInfoHsh = %{$_[0]};
	my %intronlessGeneCodingSeqHsh = %{$_[1]};
	my $minPrmntReadNumRatioAsNonStochasticAltSite = $_[2];
	my $minSenseSplicingEfficiencyAsNonStochasticIntronCreate = $_[3];

	#---get the junctStr with with complete sequence, with unq == 1 and with >0 senseSplicingEfficiency
	my $minReadNum = 1;
	my $minUnq = 1;
	my %seqForDremeAndWeblogoHsh;
	
	my $trimExonEnd = 10;
	my $trimIntronEnd = 30;
	my $lowerCutoffFactor = 10;
	
	foreach my $junctStr (keys %NGSJunctInfoHsh) {
		my $unq = ${$NGSJunctInfoHsh{$junctStr}}{"unq"};
		my $completeSeq = ${$NGSJunctInfoHsh{$junctStr}}{"completeSeq"};
		my $senseSplicingEfficiency = ${$NGSJunctInfoHsh{$junctStr}}{"senseSplicingEfficiency"};
		my $readNum = ${$NGSJunctInfoHsh{$junctStr}}{"readNum"};
		my $prmntReadNumRatio = ${$NGSJunctInfoHsh{$junctStr}}{"prmntReadNumRatio"};
		my $altSplcType = ${$NGSJunctInfoHsh{$junctStr}}{"altSplcType"};

		#----passed the filter on completeSeq , uniquness and readnum
		if (($completeSeq eq "complete") and ($unq <= $minUnq) and ($senseSplicingEfficiency ne "null") and ($readNum >= $minReadNum)) {
			
			#---determine the type of the junction
			my $junctType = 'notSelected';
			if ($altSplcType eq 'ref') {
					$junctType = 'referenceJunction';
			} elsif ($altSplcType eq 'altSite') {
				if ($prmntReadNumRatio >= $minPrmntReadNumRatioAsNonStochasticAltSite) {
					$junctType = 'nonStochasticAltSite';
				} elsif ($prmntReadNumRatio < ($minPrmntReadNumRatioAsNonStochasticAltSite/$lowerCutoffFactor)) {
					$junctType = 'stochasticAltSite';
				}
			} elsif ($altSplcType eq 'intronCreate') {
				if ($senseSplicingEfficiency >= $minSenseSplicingEfficiencyAsNonStochasticIntronCreate) {
					$junctType = 'nonStochasticIntronCreate';
				} elsif ($senseSplicingEfficiency < ($minSenseSplicingEfficiencyAsNonStochasticIntronCreate/$lowerCutoffFactor)) {
					$junctType = 'stochasticIntronCreate';
				}
			}
			
			if ($junctType ne 'notSelected') {
				my $leftFlankSeq = ${$NGSJunctInfoHsh{$junctStr}}{"leftFlankSeq"};
				my $rightFlankSeq = ${$NGSJunctInfoHsh{$junctStr}}{"rightFlankSeq"};
				my $oriLen = length $leftFlankSeq;
				my $trimLen = $oriLen-$trimExonEnd-$trimIntronEnd;
				die "the trimmed the length of the flank sequence is smaller than zero" if ($trimLen <= 0) ;

				my $trimLeftFlankSeq = substr $leftFlankSeq, $trimExonEnd, $trimLen;
				my $trimRightFlankSeq = substr $rightFlankSeq, $trimIntronEnd, $trimLen;
				my $trimFlankSeq = $trimLeftFlankSeq.$trimRightFlankSeq;

				my $leftExonSeq = substr $leftFlankSeq, 0, $boundWidth;
				my $rightExonSeq = substr $rightFlankSeq, $boundWidth+2;
				my $exonSeq = $leftExonSeq.$rightExonSeq;
				
				#my $intronSize = ${$NGSJunctInfoHsh{$junctStr}}{"intronSize"};
				my $leftIntronSeq = substr $leftFlankSeq, $boundWidth, 20;
				my $rightIntronSeq = substr $rightFlankSeq, $boundWidth - 20 + 2, 20;
				my $intronSeq = $leftIntronSeq.$rightIntronSeq;
				
				push @{${${$seqForDremeAndWeblogoHsh{$junctType}}{'intronSeq'}}{'seq'}}, $intronSeq;
				push @{${${$seqForDremeAndWeblogoHsh{$junctType}}{'exonSeq'}}{'seq'}}, $exonSeq;
				push @{${${$seqForDremeAndWeblogoHsh{$junctType}}{'trimFlankSeq'}}{'seq'}}, $trimFlankSeq;
				
				if (($junctType eq 'nonStochasticIntronCreate') or ($junctType eq 'nonStochasticAltSite')) {
					push @{${${$seqForDremeAndWeblogoHsh{'nonStochasticPool'}}{'intronSeq'}}{'seq'}}, $intronSeq;
					push @{${${$seqForDremeAndWeblogoHsh{'nonStochasticPool'}}{'exonSeq'}}{'seq'}}, $exonSeq;
					push @{${${$seqForDremeAndWeblogoHsh{'nonStochasticPool'}}{'trimFlankSeq'}}{'seq'}}, $trimFlankSeq;
				}
				
				if (($junctType eq 'stochasticIntronCreate') or ($junctType eq 'stochasticAltSite')) {
					push @{${${$seqForDremeAndWeblogoHsh{'stochasticPool'}}{'intronSeq'}}{'seq'}}, $intronSeq;
					push @{${${$seqForDremeAndWeblogoHsh{'stochasticPool'}}{'exonSeq'}}{'seq'}}, $exonSeq;
					push @{${${$seqForDremeAndWeblogoHsh{'stochasticPool'}}{'trimFlankSeq'}}{'seq'}}, $trimFlankSeq;
				}
			}
		}
	}
	
	#---print the sequences
	system ("mkdir -m 777 -p $outDir/nonStochasticDremeWebLogo/fasta/");
	system ("mkdir -m 777 -p $outDir/nonStochasticDremeWebLogo/weblogo/");
	foreach my $junctType (keys %seqForDremeAndWeblogoHsh) {
		foreach my $seqType (keys %{$seqForDremeAndWeblogoHsh{$junctType}}) {
			my $fastaPath = "$outDir/nonStochasticDremeWebLogo/fasta/$junctType.$seqType.fasta";
			my $weblogoPath = "$outDir/nonStochasticDremeWebLogo/weblogo/$junctType.$seqType.pdf";
			${${$seqForDremeAndWeblogoHsh{$junctType}}{$seqType}}{'fastaPath'} = $fastaPath;
			${${$seqForDremeAndWeblogoHsh{$junctType}}{$seqType}}{'weblogoPath'} = $weblogoPath;

			open (FASTA, ">$fastaPath");
			my $seqNum = 0;
			foreach my $seq (@{${${$seqForDremeAndWeblogoHsh{$junctType}}{$seqType}}{'seq'}}) {
				$seqNum++;
				print FASTA ">$junctType.$seqType.$seqNum\n";
				print FASTA "$seq\n";
			}
			close FASTA;
		}
	}
	
	#-----generate random sequences
	my $randCDSSeqPath = "$outDir/nonStochasticDremeWebLogo/fasta/randCDS_100nt.fasta";
	my $randCDSWeblogoPath = "$outDir/nonStochasticDremeWebLogo/weblogo/randCDS_100nt.pdf";
	my $randSeqLen = 100;
	my $randSeqNum = 1000;
	randomGenDNASeq(\%intronlessGeneCodingSeqHsh, $randSeqLen, $randSeqNum, $randCDSSeqPath);

	#-----store fake sequence in the random hash, for convenicence of retrieving
	foreach my $junctType (keys %seqForDremeAndWeblogoHsh) {
		foreach my $seqType (keys %{$seqForDremeAndWeblogoHsh{$junctType}}) {
			${${$seqForDremeAndWeblogoHsh{'random'}}{$seqType}}{'fastaPath'} = $randCDSSeqPath;
			${${$seqForDremeAndWeblogoHsh{'random'}}{$seqType}}{'weblogoPath'} = $randCDSWeblogoPath;
		}
		last;
	}

	#---define the comparison in dreme, only against neg will be made
	my %dremeComparisonPairHsh;
	$dremeComparisonPairHsh{'nonStochasticIntronCreate'} = 'stochasticIntronCreate';
	$dremeComparisonPairHsh{'nonStochasticAltSite'} = 'stochasticAltSite';
	$dremeComparisonPairHsh{'nonStochasticPool'} = 'stochasticPool';
	$dremeComparisonPairHsh{'referenceJunction'} = 'random';
	
	foreach my $querySet (keys %dremeComparisonPairHsh) {
		my $negativeSet = $dremeComparisonPairHsh{$querySet};

		foreach my $seqType (keys %{$seqForDremeAndWeblogoHsh{$querySet}}) {

			#---run weblogo for both query and negative
			foreach my $set ($querySet, $negativeSet) {
				my $fastaPath = ${${$seqForDremeAndWeblogoHsh{$set}}{$seqType}}{'fastaPath'};
				my $weblogoPath = ${${$seqForDremeAndWeblogoHsh{$set}}{$seqType}}{'weblogoPath'};
				my $numSeq = `sed -n \'\$=\' $fastaPath`;
				$numSeq = int ($numSeq/2);
				my $title = "$set.$seqType";
				$title .= "_n=$numSeq";
				my $length = 150;

				print "Running weblogo for $set $seqType\n";
				my $weblogoCMD = "weblogo -f $fastaPath -F pdf -n $length -s large --title $title -A rna -c classic >$weblogoPath";
				my $grepCMD = "ps -ef | grep weblogo | grep -v grep | grep -v perl";
				runAndCheckSerialTask($grepCMD, "weblogo", $weblogoCMD, "$outDir/error.log.txt");
			}
		
			#---run dreme
			my $queryFastaPath = ${${$seqForDremeAndWeblogoHsh{$querySet}}{$seqType}}{'fastaPath'};
			my $negativeFastaPath = ${${$seqForDremeAndWeblogoHsh{$negativeSet}}{$seqType}}{'fastaPath'};
		
			my $dirForDremePath = "$outDir/nonStochasticDremeWebLogo/dreme/$seqType.$querySet.vs.$negativeSet/";
			system "mkdir -p -m 777 $dirForDremePath";
			my $dremeCMD = "dreme -oc $dirForDremePath -p $queryFastaPath -n $negativeFastaPath -maxk 10 -norc -mink 4 -e 0.001";
			my $grepCMD = "ps -ef | grep dreme | grep -v grep | grep -v perl";
			my $grepStr = "dreme -oc $dirForDremePath";
			runAndCheckSerialTask($grepCMD, $grepStr, $dremeCMD, "$outDir/error.log.txt");
			plotDremeMotif($dirForDremePath, $dirForDremePath."./dreme.xml", $queryFastaPath, "$querySet.$seqType", "$seqType of $querySet vs $negativeSet");
		}
	}
}
