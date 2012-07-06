#########################################################################
######	TO DO LIST:
#		i	make summary compare algorithm
#		ii	build package with package.skeleton
#		iii	
#########################################################################
#########################################################################
######	UPDATE LOG - please add new contributions in the format:
######
#########################################################################
#########################################################################
#[date]	i	what you did
#		ii	second thing you did
#		iii	etc
#########################################################################
#########################################################################
#062112		i 	made cutsites for msgfdb output  much quicker 
#
#
#061512	i	fixed list.peps, unique.peps, and unique.peps so they are faster
#
#061312	i 	building a function to plot success rates of identificatoins as a function of m/z of precursor
#			this function requires the .mgf input file and the .txt MSGFDB output file
#
#061212	i	broke up pepsum into multiple functions, created "calibrate.mgf()" 
#			which is more flexible than the previous method and much faster
#
#042412	i	started making a pepsum function to summarize and update the spectra
#		ii 	pepsum function can update the precursor and fragment masses if the 
#			spectra are in .mgf format
#
#041312	i	added function normalizetable() to divide each position in a table by the sum of the table
#
#
#041212	i 	added function to perform peptide string parsing of MSGFDB peptide strings
#		ii	added @pepvec and @provec to pepsum object slots
#
#040512	i	read.MSGFDB input fixed such that
#			now the total and the total-by-fraction are printed
#
#032812	i	fixed cut.predict[works with n-term or c-term] and coverage 
#
#031812	i	finished the functions spectra.count and genQspec to automatically 
#			generate spectral counts and a Qspec input file for spectral counting
#			the file should be uploaded to http://www.nesvilab.org/qspec.php/
#			the genQSpec output file needs to have the ""'s removed manually before upload
#
#122111	i	indexfiles now allows indexing by activation, could easily add others
#			format of index is list with named slots containing index numbers allowing legacy support for my previous indexing
#		ii	read.MSGFDB outputs an update about which file it was on
#		iii	new code base: get.peps added:  uses new indexing to get specific subsets meeting multiple criteria from  MSGFDB objects
#
#122011	i	updated indexfiles to allow indexing by protein accession, fraction, and pval for MSGFDB type files
#
#
#121311	i	added companion functions to calculate monoisotopic mass
#			and expand those M values into m/z for M+1H to M+5H
#		ii	nearly complete with the area compare interface through xcmsEIC class usage
#
#120311	i	fixed read.proteome to deal with zeros
#		
#		i	decided to combind the proteome into the pepsum object
#			to allow for normalized spectral counting
#			WORKS but not a function yet..
#
#110711	i	started adding a function to parse data by fraction and confidence to
#			output a _syn.txt type annotation file for usage with SpectrumLook 
#			program from PNNL that shows annotations given .mzXML and .txt input
#
#			
#
#092111	i	added the function spectra.count and it works for MSGFDB type objects
#			need to change this so that it makes a csv instead with the type:
#			For example (let's call it "mydata.csv"), sets C and W, each with 4 replicates:
#			IPI,C4,C3,C2,C1,W4,W3,W2,W1
#			IPI00323624,106,118.5,103,65,64.5,67.5,52.5,49.5
#			IPI00553798,78,57.5,81.5,69,67,67,64.5,94.5
#			IPI00131138,82,87.5,87,87,89.5,77,91,85.5
#			IPI00227299,87.5,80,90,95.5,158,100.5,153,107
#			
#090811	i	added @summary to pepsum, read.prospector runs Pep.suM before exit
#
#090411	i	Pep.suM now works with matrix output of search stats
#		ii	Pep.meR started to compare significant peptides across search engines
#		iii	multiple other function improvements for flexibility
#		iv 	charge count works!!! plots charge vs. basic residue 
#			with linear correlation analysis
#		v	started EIC.calc
#		
#072511	i	Added read.MSGFDB to functions and it works great
#
#
#061011	i	cutsites now outputs .jpeg with cutsites(title=""), allows multiple 
#			pictures to be generated
#
#
#051511	i 	started class PROsum for protein summary files with iTRAQ quan
#
#051411	i 	started making a function to spit out protein names for 
#			object@filetype="pilot"
#		ii	structured this file so anyone can work on it
#
#050911	i	read.prospector(index=TRUE) now works and makes 
#			object@index=matrix where columns are indexes of line #s 
#			which correlate to each fraction in the .txt source
#			This will allow the rapid parsing of fraction lines 
#			while generating cutsites heatmaps 
#
#
#
#050811	i	make object@index into a list with indexes correlated to name of the list entry; e.g. fraction1-s1<-c()
#
#050711	i	cutsites(exclude=c()) now works with character-type fractions and "prospector" type pepsum objects 
#
#050411	i	read.pilot now works with multiple files or a single file fed as "input" 
#			either input generates a data.frame with the file name in the last row
#		ii	read.inspect works and cutsites can now process data from prospector, pilot and inspect!!!!
#
#
#
#040311	i	fixed conditional prospector versus pilot in cutsites()
#		ii	cutsites() with prospector is lightning fast, truncated steps
#
#033011	i	started the indexfiles()  for prospector
#
#032911	i	fixed read.prospecter to work with the union of searches merged
#
#032311	i	added genomestats()
#		ii	fixed (i hope) cutsites
#
#
#013011
#		i	started read.mascot
#		ii	started read.ommsa
#
#
#012911
#		i	separated the filtering and indexing into distinct functions that work on
#			objects of class pepsum, works for read.prospector 
#		ii	made filterfiles work with object@filetype="pilot" checks is.null(object@index)
#			loops through indexed numbers and builds object@sequence 
#		iii	cutsites now works with objects of filetype "pilot" and "prospector"
#
#012711	1	read.pilot works and indexes peptide lines that are above conf
#
#
#########################################################################
####	#	#	#	#	#	#	#	#	#	#	#	######
####	CLASSES AND ENVIRONMENTS FIRST						######
####	define class pepsum, make blank object of class pepsum		######
#### 	define amalgamation of pepsum and xcmxSet of raw data 		######
####	files for co-manipulation							######
####	#	#	#	#	#	#	#	#	#	#	#	######
#########################################################################
#########################################################################
####	New environment for objects to be accessible				######
####												######
#########################################################################

PEPSUM<-new.env(hash = FALSE, parent = parent.frame(),obj1)
require(XML)
require(xcms)

#setClass("PepMain",representation(pepsum="pepsum",xcmsSet="xcmsSet",proteome="proteome"))

setClass("pepsum",representation(summary="matrix",scans="ANY",specdir="character",specfiles="list",data="ANY", totalheat="matrix",residues="character",
	sequence="ANY",score="character",filetype="character", index="ANY", fraction="character",pepvec="ANY",provec="ANY",filenames="ANY",p1.index="ANY",countlist="list"))


setClass("proteome",representation(fasta="list", AAratios="matrix",seq="list",names="ANY",cleaved="list",coverage="ANY",summary="ANY"))


object<-new("pepsum")

setClass("prosum",representation(data="ANY",filetype="character",index="ANY",fraction="character",filenames="ANY",proteins="character"))



###Function leader format below for when the package is ready to build
#########################################################################
######	Does it work?   yes/no!!!!! [date]
#########################################################################
######	Function does: what????
######	useful when:  when?????
#########################################################################
#####		usage:  example usage???
#########################################################################




#########################################################################
#########################################################################
#################
####function to read fasta file and make object with stats
#########################################################################


read.proteome=function(input=fastapicker(),dir=getwd(),cleave=TRUE,pie=FALSE){
    setwd(dir)
    ##### read .fasta file and fill up object
    #####    slots with useful tidbits
    require(seqinr)
    require(gplots)
    object<-new("proteome")
    object@fasta<-read.fasta(input,seqtype="AA",as.string=TRUE)
    object@seq<-getSequence(object@fasta)
    object@names<-names(object@fasta)
    x<-unlist(object@seq)
    aatotal<-length(x)
    ###empty summary matrix
    summary<-matrix(c(rep(0,40)),nrow=20,dimnames=list(
        c("A","S","T","G","V","C","N","L","I","M","P","Y","W","Q","F","D","E","H","K","R"),
        c("count","ratio")),
        ncol=2)
	summary[,1]<-c(length(x[x=="A"]),length(x[x=="S"]),length(x[x=="T"]),
      length(x[x=="G"]),length(x[x=="V"]),length(x[x=="C"]),
	length(x[x=="N"]),length(x[x=="L"]),length(x[x=="I"]),
	length(x[x=="M"]),length(x[x=="P"]),length(x[x=="Y"]),
	length(x[x=="W"]),length(x[x=="Q"]),length(x[x=="F"]),
	length(x[x=="D"]),length(x[x=="E"]),length(x[x=="H"]),
	length(x[x=="K"]),length(x[x=="R"]))
	for(i in 1:length(summary[,1])){summary[i,2]<-summary[i,1]/aatotal}
	par(mfcol=c(1,2))
	if(pie==TRUE){
		pie(summary[,2],main=paste(input,"stats"),radius=1)
		}
	textplot(summary)
	object@AAratios<-summary
	print(aatotal)
	if(cleave==TRUE){
		object<-cut.predict(object=object,residues=c("R","K"))
		}
	return(object)
	} 
 

#############################################################################
######		use the above class "proteome" to predict cutsites
######		given a single amino acid
#############################################################################
######	082811	 this works
#############################################################################

################################################################################
## cut predictor v2; for use with multiple cleavage points
###########################
####				9/1/2011 THIS WORKS with an arbitrary number of residues

##### version three, allow cleavage at n-terminal or c-terminal
####	use a vector of length == length(residues) where
####	0 equals n-terminal to the residue, and 1 equals c-terminal to the residue
#####	
####### 032712 this works to cleave at either n-term or c-term
##########		#############		###########
###	FIXED: the loop when term=0 needs to start on 2 instead of 1


cut.predict=function(object=spombe,residues=c("D"),term=c(0),protease="AspN"){
	###	testing the cleavage n-terminal to residues[1]
	#########
	proteins.l<-length(object@fasta)
	peps<-c(rep(0,sum(nchar(unlist(object@seq)))))
	peps2<-c(rep(0,sum(nchar(unlist(object@seq)))))
	pep.lengths<-c(rep(0,sum(nchar(unlist(object@seq)))))
	peps.l<-length(peps)
	lengths<-c()
	print(residues[1])
	m<-1###the location in peps
	for(i in 1:proteins.l){
		pro<-unlist(object@seq[i])
		pro.l<-length(pro)
		###loop over AAs in pro[i]
		k<-1
		if(term[1]==0){
			for(j in 2:pro.l){
				if(pro[j]==residues[1]){
					end<-j-1
					peps[m]<-c(paste(pro[k:end],collapse=""))
					k<-j
					m<-m+1
					}
				if(j==pro.l){
					peps[m]<-c(paste(pro[k:j],collapse=""))
					m<-m+1
					}
				}
			}
		##### if cleavage is supposed to be at the n-terminal 
		if(term[1]==1){
			for(j in 1:pro.l){
				if(pro[j]==residues[1]){
					peps[m]<-c(paste(pro[k:j],collapse=""))
					k<-j+1
					m<-m+1
					}
				else if(j==pro.l){
					peps[m]<-c(paste(pro[k:j],collapse=""))
					m<-m+1
					}
				}
			}
		}
	if(length(residues)>1){
		for(n in 2:length(residues)){
			print(residues[1:n])
			peps.l<-length(peps[peps!=0])
			m<-1
			for(i in 1:peps.l){
        			pep<-unlist(strsplit(peps[i],split=""))
				pep.l<-length(pep)
				###loop over AAs in pro[i]
				k<-1
				for(j in 1:pep.l){
					if(pep[j]==residues[n]){
					##if the jth residue=residues
					###write pep 1:j
						peps2[m]<-c(paste(pep[k:j],collapse=""))
						k<-j+1
						m<-m+1
						}
					###if the jth residue == the last residue
					else if(j==pep.l){
						peps2[m]<-c(paste(pep[k:j],collapse=""))
						k<-j+1
						m<-m+1
						}
					}
				}
			peps<-peps2
			}
		}
	peps<-peps[peps!=0]
	object@cleaved<-list(protease<-peps)
	return(object)
	}


######################################################################################
#################################################################
########		script to find coverage given cleavage output from above
##########################################################################
#########	adapt this to combine length(object@cleaved)

coverage=function(object=Spombe){
	#### how many proteases were used?
	cleaved.l<-length(object@cleaved)
	###loop over the cleaved.l
	for(j in 1:cleaved.l){
		peps2<-unlist(object@cleaved[j])
		totalvec<-c(rep(0,sum(nchar(peps2))))
		peplengths<-nchar(peps2)
		pepvec.l<-length(peplengths)
		start<-1
		end<-peplengths[1]
		if(peplengths[1]>=7&peplengths[1]<=35){
			totalvec[start:end]<-rep(1,times=peplengths[1])
			}
		start<-end+1
		for(i in 2:pepvec.l){
			end<-peplengths[i]+end
			if(peplengths[i]>=7&&peplengths[i]<=35){
				totalvec[start:end]<-rep(1,times=peplengths[i])
				}
			start<-start+peplengths[i]
			print(i)
			}
		### this part puts the binary outcome into object@coverage
		if(j==1){object@coverage<-totalvec}
		if(j>=2){object@coverage<-object@coverage+totalvec}
		}
	notcovered<-length(totalvec[totalvec==0])
	total<-length(totalvec)
	pctnotcovered<-notcovered/total
	pctcovered<-1-pctnotcovered
	object@summary<-c(covered=pctcovered,not=pctnotcovered)
	pie(summary,main="predicted coverage")
	return(object)
	}



#########################################################################
######	read.MSGFDB works and prints unique peptides and protein counts
######	can read a directory of input files done
#########################################################################
######	Function to generate object of class pepsum from MSGFDB
######	output from tabdelim.txt files
######	index=TRUE indexes row numbers for <=0.01 Pepfdr 
#########################################################################
######	usage: read.MSGFDB(dir=getwd(),input="file.txt",index=TRUE)
######	unique=TRUE prints # of unique peptides and proteins, respectively
#########################################################################

read.MSGFDB=function(dir=getwd(),input=list.files(),index=TRUE, unique=TRUE,by=list(pval=0.01),pepsum=TRUE){
	object<-new("pepsum")
	object@filetype="MSGFDB"
	object@data<-read.delim(input[1],fill=TRUE,header=TRUE,row.names=NULL)
	object@fraction<-input
	source.col<-length(object@data)+1
	object@data[,source.col]<-input[1]
	print(1)
	if(length(object@fraction)>=2){
		for(i in 2:length(object@fraction)){
			u<-read.delim(input[i],header=TRUE)
			source.col<-length(u)+1
			u[,source.col]<-input[i]
			object@data<-merge(u,object@data,all=TRUE,sort=FALSE)
			print(i)
			}
		}
	if(index=="TRUE"){object<-indexfiles(object,by=by)}
	if(unique=="TRUE"){
		print("total unique peptides")
		object<-unique.peps(object)
		print(length(object@pepvec))
		print("unique peptides by fraction")
		uniquepeptidesbyfrac(object)
		print("unique Proteins by fraction")
		unique.proteinsbyfrac(object)
		print("total unique proteins")
		object<-unique.pros(object)
		print(length(object@provec))
		print("# of PSMs")
		print(length(list.peptides(object)))
		}

	if(pepsum){pepsum(object)}
	return(object)
	}

#########################################################################
######		addition of fragmentation method to params works
######	unique.peps AND unique.pros for use with read.MSGFDB()
#########################################################################
######	Function to generate list of unique peps from MSGFDB output
######	object@filetype==MSGFDB; object@index$pepfdrs
######		only works by all fractions
#########################################################################
#####		usage: unique.peps(object)
#########################################################################


unique.peps=function(object){
	if(object@filetype!="MSGFDB"){
		print("wrong file type")
		}
	if(object@filetype=="MSGFDB"){
		sigpeps<-object@index$pepfdrs[object@index$pepfdrs!=0]
		peps<-as.character(object@data$Peptide[sigpeps])
		peps<-unique(peps)
		peps<-na.omit(peps)
		object@pepvec<-peps
		}
	return(object)
	}


unique.pros<-function(object){
	if(object@filetype!="MSGFDB"){
		print("wrong file type!")
		}	
	if(object@filetype=="MSGFDB"){
		sigpeps<-object@index$pepfdrs[object@index$pepfdrs!=0]
		pros<-as.character(object@data$Protein[sigpeps])
		pros<-unique(pros)
		pros<-na.omit(pros)
		object@provec<-pros
		}
	return(object)
	}

unique.proteinsbyfrac=function(object=alp){
	for(i in 1:length(object@index$fractions[1,])){
		fracPval<-object@index$fractions[,i]+object@index$binaryfdr
		uniqueProteins<-length(unique(object@data$Protein[grep(fracPval,pattern=2)]))
		print(object@fraction[i])
		print(uniqueProteins)
		}
	}

#########################################################################
######	040512	
######	for use with read.MSGFDB()
#########################################################################
######	Function to generate list of unique peps from MSGFDB output
######	object@filetype==MSGFDB; object@index$pepfdrs
######		only works by all fractions
#########################################################################
#####		usage: unique.peps(object)
#########################################################################

uniquepeptidesbyfrac=function(object=alp){
	for(i in 1:length(object@index$fractions[1,])){
		fracPval<-object@index$fractions[,i]+object@index$binaryfdr
		uniquePeptides<-length(unique(object@data$Peptide[grep(fracPval,pattern=2)]))
		print(object@fraction[i])
		print(uniquePeptides)
		}
	}


unique.pepsbyfrag<-function(object,FragMethod="ETD"){
	if(object@filetype=="MSGFDB"){
		peps<-c()
		pros<-c()
		i<-1
		for(x in object@index$pepfdrs){
			if(x!=0){
			#print(as.character(object@data[x,"FragMethod"]))
				if(as.character(object@data[x,"FragMethod"])==FragMethod){
					peps[i]<-as.character(object@data$Peptide[x])
					pros[i]<-as.character(object@data$Protein[x])
					i=i+1
					}
				}
			}
		print("PSMs")
		print(length(peps))
		#peps<-unique(peps)
		#pros<-unique(pros)
		print("unique peptide count")
		print(length(unique(peps)))
		print("unique protein count")
		print(length(unique(pros)))

		both<-matrix(c(peps,pros),ncol=2,byrow=FALSE)
		return(both)
		}
	if(object@filetype!="MSGFDB"){
		print("wrong file type")
		}
	}


######################################################################################################
#########################################################################
######	092811	the below functions allow protein mod finding
######			similar to the above except not unique only
#########################################################################
######	NOT finished
######	
######	
#########################################################################
#####		usage: unique.peps(object)
#########################################################################
list.proteins<-function(input=object){
	object=input
	pros<-c()
	if(object@filetype=="MSGFDB"){
		for(x in object@index$pepfdrs)pros[x]<-as.character(object@data$Protein[x])
		}
	if(object@filetype!="MSGFDB"){
		for(x in object@index$pepfdrs)pros[x]<-as.character(object@data$Protein.Name[x])
		}
	pros<-na.omit(pros)
	return(pros)
	}

list.peptides<-function(input=object){
	object=input
	pepfdrs<-object@index$pepfdrs
	sigpeps<-pepfdrs[pepfdrs!=0]
	if(object@filetype=="MSGFDB"){
		peps<-as.character(object@data$Peptide[sigpeps])
		}
	if(object@filetype!="MSGFDB"){
		for(x in object@index$pepfdrs)peps[x]<-as.character(object@data$DB.Peptide[x])
		}
	peps<-na.omit(peps)
	return(peps)
	}





list.PMError<-function(input=object){
	object=input
	pme<-c()
	if(object@filetype=="MSGFDB"){
		for(x in object@index$pepfdrs)pme[x]<-as.character(object@data$PMError.ppm.[x])
		}


	if(object@filetype=="prospector"){
		for(x in object@index$pepfdrs)pme[x]<-as.character(object@data$ppm[x])
		}
	
	pme<-na.omit(pme)
	return(pme)
	}


list.lines<-function(object){
	if(object@filetype=="MSGFDB"){
		sigpeps<-object@index$pepfdrs
		lines<-object@data[sigpeps,]
		}
	if(object@filetype=="prospector"){
		print("prospector doesnt work yet")
		}
	lines<-na.omit(lines)
	return(lines)
	}

findmods2<-function(object,modmass="79",write=FALSE,output="sumos.txt"){
	lines<-list.lines(object)

	indexed<-grep(lines[["Peptide"]],pattern=modmass)###gives the indexes of the peptides with +114 mods
	lines<-lines[indexed,]

	return(lines)
	if(write){
		write.table(lines,file=output,sep="\t", row.names = F, col.names = F)
		}
	}

##############################################################
####		this function will take a pepsum object of type MSGFDB
###		summarizes the error
####		use this to update the spectral precursor masses
#####	04/24/2012 	This appears to work and outputs a new mgf
####		with modified precursor masses
###			also updates the fragment ions 
#### the calibration only works for .mgf spectral files
#######################################

pepsum=function(object,modmass=c(114),cal=TRUE,specfile=c("C:/MSGFDB.20120106/data/sumo/tS.mgf"),name="aLP_CID"){
	sigpeps<-object@index$pepfdr
	peps<-unique(as.character(object@data[sigpeps,"Peptide"]))
	masserror<-object@data[sigpeps,"PMError.ppm."]
	meanerror<-sum(masserror)/length(masserror)
	hist(masserror,breaks=seq(min(masserror)-5,max(masserror)+5,by=0.1),main=name)
	abline(v=meanerror,col="blue",lwd=2)
	if(length(modmass)>=1){
		mods<-findmods(object)
		print(mods)
		}
	
	###text(labels="mass error"
	#### this part takes the spec file and updates the precursor masses
	### there is a problem due to the \t between mass and intensity
	if(cal & length(specfile)>=1){
		calibrate.mgf(object, specfile=specfile, name=name)
		}
	}


#####  	his function will take a pepsum object and its mgf spectra file 
#####		and updates the masses using the observed mean 
#####		outputs a new file with "cal" appended to it
#####		this works but it is not flexible



#### almost works
#####
split.mgf=function(object,specfile=c("alp_CIDHCDETDSA.mgf"),name="aLP_CIDHCDETDSA",all=TRUE,activation=""){
	if(all){
		con<-file(description = specfile, blocking = FALSE, open = "r")
		alllines<-readLines(con)
		close(con)
		activationlines<-grep("ACTIVATION=",alllines)
		activationlines.l<-length(activationlines)
		activations<-alllines[activationlines]
		activationlevs<-levels(as.factor(activations))
		activationlevs.l<-length(activationlevs)
		lengthperactivation<-c()
		for(x in activationlevs){
			lengthperactivation<-c(lengthperactivation,length(activations[activations==x]))
			}
		for(x in activations){
		##### 	make new files with length equal to each 
		######	and then fill up those lines into the new spectra



		}


	if(all==FALSE){
		}
	}

calibrate.mgf=function(object=alp5ppmCIDETD, specfile=c("C:/msgfdb/data/alp_velos_HPRP.mgf"),name="aLP_velos_HPRP"){
	sigpeps<-object@index$pepfdr
	masserror<-object@data[sigpeps,"PMError.ppm."]
	meanerror<-sum(masserror)/length(masserror)
	hist(masserror,breaks=seq(min(masserror)-5,max(masserror)+5,by=0.1),main=name)
	abline(v=meanerror,col="blue",lwd=2)
	con<-file(description = specfile, blocking = FALSE, open = "r")
	alllines<-readLines(con)
	close(con)
	masslines<-grep("PEPMASS=",alllines)
	masslines.l<-length(masslines)
	mass<-as.numeric(substr(alllines[masslines],9,17))
	### make a new vector of updated masses
	newmass<-mass-mass*(meanerror/1000000)
	newmass<-round(newmass,digits=5)
	newlines<-paste("PEPMASS=",newmass,sep="")
	newfile<-alllines
	newfile[masslines]<-newlines
	###  this part is for updating product ion masses
	### index lines where mass, intensity pairs start and end
	startions<-grep("BEGIN IONS",alllines)+6
	endions<-grep("END IONS",alllines)-1
	spectra.l<-length(startions)
	### the positions need to change based on the file lines
	### in some cases 5 , 6 or 7 lines
	####positions<-rep(0,length(alllines)-length(startions)*7)
	for(i in 1:spectra.l){
		if(startions[i]<endions[i]){
		posvec<-seq(startions[i],endions[i],by=1)
		posvec.l<-length(posvec)
		lines<-alllines[posvec]
		mzpos<-seq(1,posvec.l*2,by=2)
		intpos<-seq(2,posvec.l*2,by=2)
		mzint<-unlist(strsplit(lines,split="\t"))
	### as.numeric works but the numbers get rounded to 4 decimals
		mz<-as.numeric(mzint[mzpos])
		mzint[mzpos]<-as.character(round(mz-mz*(meanerror/1000000),digits=4))
		newfile[posvec]<-paste(mzint[mzpos],mzint[intpos],sep=" ")
		print(i)
		}
		}
	length(newfile)
	length(alllines)
	cat(newfile, file = paste(name,"cal.mgf",sep="_"), sep = "\n", fill = FALSE, labels = NULL,
		append = FALSE)
	}
	
calibrate.mgf()

##################################################################################################
########################    this works but not as a function 		##########################
###############################		06/14/12			################################
##########################							    ############################
##################################################################################################

plotmzefficiency=function(object=alp5ppmCIDETD, chunks=20,pval=0.001,specfile=c("C:/msgfdb/data/alp_cidhcdetdsa_cal.mgf"),name="aLP_CIDETD"){
	#### how many activations
	activationlevs<-levels(as.factor(object@data[,"FragMethod"]))
	activationlev.l<-length(activationlevs)
	fragvector<-object@data[,"FragMethod"]
	activationlengths<-c()
	#### how many per activation
	for(i in 1:activations.l){
		activationlengths<-c(activationlengths,length(fragvector[fragvector==activations[i]]))
		}
	### how many per charge at each activation?
	zlevs<-as.numeric(levels(as.factor(object@data[,"Charge"])))
	z.l<-length(zlevs)
	#### get the output specific range of data
	precrange<-range(object@data[,"Precursor"])
	#### divide the range into 10 segments
	##### this part is to get the identifications from the MSGFDB output
	##### take all columns from the significant rows 
	sigpeps<-object@index$pepfdr
	identified<-object@data[sigpeps,]
	IDscans<-identified[,"Scan."]
	#### might not need to go into the mgf actually because the output has all spectra
	#### part to get the precursor info from the .mgf file
	con<-file(description = specfile, blocking = FALSE, open = "r")
	alllines<-readLines(con)
	close(con)
	activationlines<-grep("ACTIVATION=",alllines)
	activations<-alllines[activationlines]
	actchar<-substr(activations,start=12,stop=14)
	masslines<-grep("PEPMASS=",alllines)
	masses<-alllines[masslines]
	massints<-as.numeric(substr(masses,start=9,stop=17))
	scanlines<-grep("SCAN",alllines)	
	scans<-alllines[scanlines]
	scanintegers<-as.numeric(substr(scans,start=7,stop=11))
	chargelines<-grep("CHARGE",alllines)
	charges<-alllines[chargelines]
	chargeints<-as.numeric(substr(charges,start=8,stop=8))
	### get rid of scans where the charge was undetermined
	ids<-rep(0,times=length(scans))
	#### loop through the identified scans and put 1's into those scan positions
	### works but kinda slowly	
	IDscanlen<-length(IDscans)
	i=0
	for(x in IDscans){
		ids[grep(scanintegers,pattern=x)]=1
		i=i+1
		print(i/IDscanlen)
		}
	idframe<-data.frame(scanintegers,actchar,chargeints,massints,ids)
	idorder<-idframe[ order(idframe[,4]),]
	spectra.len<-length(idorder[,1])
	#### need empty lists for each
	## divide the ordered matricies into their charges
	bycharge<-list()
	for(i in 1:z.l){
		bycharge[[i]]<-idorder[grep(idorder[,3],pattern=zlevs[i],fixed=T),]
		}
	### divide those by activation
	byactivation<-list()
	for(i in 1:activationlev.l){
		byactivation[[i]]<-bycharge[[3]][grep(as.character(bycharge[[3]][,2]),pattern=activationlevs[i]),]
		}
	#### now that charge value is split into activations
	#### compute success as a function of bin
	#### this part works to plot the success rates 
	chunks<-15
	chunksize<-(precrange[2]-precrange[1])/chunks
	x<-seq(from=precrange[1]+chunksize,to=precrange[2], by= chunksize)
	#### which tempframe are we using?
	tempframe<-byactivation[[3]]
	mzlim<-precrange[1]+chunksize
	pcts<-rep(0, times=chunks)
	j<-0
	i=1
	for(k in 1:chunks){
		sum=0
		while(tempframe[i,4]<= mzlim & i<=length(tempframe[,4])){
			sum=sum+tempframe[i,5]
			i=i+1
			}
		pcts[k]<-sum/(i-j)
		j<-i
		mzlim<-mzlim+chunksize
		}
	hcd2<-pcts
	etd2<-pcts
	cid2<-pcts
	plot(x, hcd2,col="black",xlab="m/z",ylab="ratio of identifications",pch=3)
	lines(x,hcd2,col="black")
	title("15 chunks, +4")
	points(x,etd2,pch=2,col="red")
	lines(x,etd2,col="red")
	points(x,cid2,pch=2,col="blue")
	lines(x,cid2,col="blue")
	legend(x=450,y=0.15,legend=c("hcd","etd","cid"),pch=c(3,2,1),col=c("black","red","blue"))






##############################
###this function appears to work and outputs a tab delim .txt in descending order
####works but not for making several with fractions

####this may actually work for multiple fraction pepsum objects?

spectra.count<-function(input=object,table=TRUE){	
	object=input
	### only works for msgfdb type files
	if(object@filetype=="MSGFDB"){
		###loop through the significant PSMs in each fraction	
		fractions.l<-length(object@index$fractions[1,])
		countlist<-list()
		for(j in 1:fractions.l){
			fractionindex<-object@index$fractions[,j]+object@index$binaryfdr
			fractionindex<-grep(fractionindex, pattern=2)

			###this is super innefficient!!!!!
			######change this to use predicted empty vector memory blocks
			######empty vectors
			pros<-c()
			fragmethod<-c()

			#####	loop through the index of peptides below pval in each fraction
			for(x in fractionindex){
				pros[x]<-as.character(object@data$Protein[x])
				}
			#pros<-unique(pros)
			pros<-na.omit(pros)		
			counts<-c()
			unique.proteins<-unique(pros)
			proteins.l<-length(unique.proteins)

			###index the locations of the unique proteins
			#proindex<-c()
	
			#for(i in 1:proteins.l){
			#	proindex[i]<-length(pros[pros==unique.proteins[i]])
			#	}

			for(i in 1:proteins.l){
				counts[i]<-length(pros[pros==unique.proteins[i]])
				}
			x<-data.frame(proteins=unique.proteins,counts=counts)
			ordered<-x[order(x[,2],x[,1],decreasing=TRUE),]
			print(j)
			### keep ordered in an object
			countlist[j]<-list(ordered)
			}
		}
	if(table==TRUE){
		write.table(file="ordered.txt",sep="\t",ordered)
		}
		if(object@filetype!="MSGFDB"){
		print("wrong file type!")
		}

	###make a slot in pepsum for the spectral counts
	#### name each position in the list by its file
	names(countlist)<-object@fraction
	####output the list of counts
	object@countlist<-countlist
	return(object)
	}

#mouse<-read.proteome(input="c:/msgfdb/database/mouse.cc.fasta",cleave=FALSE)

######################################################################################################################################################################
###################################################################################
#####################[][][][][][][]################################################
###################################################################################


#### need to make a different function based on what Qspec wants as input

###prot ID		Prot Len		0		0		1	1	
###make it accept pepsum objects
###control and treatment tell which to assign ones and zeros to
#######################################################################################
#############	this works but not automatically yet...


##### control and treatment are equal to the positions in the msgfdb object
####################################
#########################################################################
#### this should work but has not been tested
####fixed this on 3/26/2012
################################

genQSpec=function(object,treatment=c(1,2,7,8),control=c(3,4,5,6),proteome=mouse,title="cavOEagedVSyoung"){
	####if countlist is empty, then run spectra.count on the object##
	if(length(object@countlist)<=1){
		object<-spectra.count(input=object,table=FALSE)
		}
	####else just run the gen
	countlist<-object@countlist
	both<-c(treatment,control)
	both.l<-length(both)

	####start with the first two
	### merged adds the x value to the left of the y value columns
	aligned<-merge(x=countlist[[both[1]]],y=countlist[[both[2]]],all=T,by=1)
	for(i in 3:both.l){
		aligned<-merge(y=aligned,x=countlist[[both[i]]],all=T,by=1)
		}

	## replace NAs with 0
	aligned[is.na(aligned)]<-0
	### dig through the proteome object for the hits
	hits<-as.character(aligned[,1])
	hits.l<-length(hits)
	index<-rep(0,hits.l)
	for(i in 1:hits.l){###NOTE: must use fixed equals true otherwise it returns everything
		index[i]<-grep(pattern=hits[i],proteome@names,value=FALSE,fixed=TRUE)
		}
	lengths<-rep(0,hits.l)
	for(j in 1:hits.l){
		lengths[j]<-length(unlist(strsplit(unlist(proteome@fasta[index[j]]),split="")))
		}
	####loop over the input to build the final .tsv file
	final<-cbind(aligned[1],lengths)
	both.l<-length(both)
	for(i in 1:both.l){
		i=i+1
		final<-cbind(final,aligned[i])
		}
	### the primitive way in case the above fails
	##final<-cbind(allaligned[1],lengths,allaligned[2],allaligned[3],allaligned[4],allaligned[5],allaligned[6],allaligned[7],allaligned[8],allaligned[9])
	colnames(final)<-c("protid","protLen",rep(0,times=length(control)),rep(1,times=length(treatment)))
	write.table(x=final,file=paste(title,"Qspec.tsv"),sep="\t",row.names=F)
	return(final)
	}


####currently this function takes a single count data frame and proteome 
#####  and normalizes based on the lengths found in the protein

###can I use this function to tap[ply my data transform?????
################################
############
##################

norml.ize=function(counts=m142n1counts,proteome=mouse){	
	###the counts data frame, values to find
	hits<-as.character(counts[,1])
	###make the proteome object with the sequences tofind
	
	hits.l<-length(hits)
	index<-rep(0,hits.l)
	for(i in 1:hits.l){
		###NOTE: must use fixed equals true otherwise it returns everything
		index[i]<-grep(pattern=hits[i],proteome@names,value=FALSE,fixed=TRUE)
		}
	lengths<-rep(0,hits.l)
	for(j in 1:hits.l){
		lengths[j]<-length(unlist(proteome@fasta[index[j]]))
		}

	spectralabundance=counts[,2]/lengths
	normsum<-sum(spectralabundance)
	NSAF=spectralabundance/normsum
	par(mfcol=c(1,2))
	hist(NSAF,breaks=1000,main="noramized spectral abundance factor")
	hist(log(NSAF),breaks=500,main="log transformation")

	
	sp.c<-data.frame(counts[1],counts[2],NSAF)
	return(sp.c)

	}



###################################################################################
####mouse###	a function to submit multiple pepsum objects to the above function 
#######	and generate a pepC compatible .csv for pepc.jar
#######	#######################################################################
###################################################################################
pepc.generator=function(input<-list(n72.ordered,n73,wt1.ordered,wt3.ordered))
	all.proteins<-c()
	input.l<-length(input)
	for(i in 1:input.l){
		all.proteins<-c(all.proteins,as.character(input[[i]]$proteins))
		}
	unique.proteins<-unique(all.proteins)
	pepc.mat<-matrix(ncol=length(input),nrow=length(unique.proteins))
	for(i in 1:input.l){
		dataframe<-unlist(input[i])
		for(j in 1:length(dataframe){
			pepc.mat[x,i]<-
			}
		}
	dataframe<-input[1]	
	typeof(dataframe[[1]])
	}

##########below can do only ETD or CID if you spec.
unique.pros2<-function(input=object,FragMethod="ETD"){
	object=input
	if(object@filetype=="MSGFDB"){
		pros<-c()
		for(x in object@index$pepfdrs)if(object@data$FragMethod[x]==FragMethod)pros[x]<-as.character(object@data$Protein[x])
		pros<-unique(pros)
		return(pros)
	}
	if(object@filetype!="MSGFDB"){
		print("wrong file type!")
	}
}


#########################################################################
####	013011 read.Inspect								#
####												#
#########################################################################


read.inspect=function(dir=getwd(),input=list.files(),index="FALSE"){
	require(gdata)
	setwd(dir)
	####make new object and start filling it################
	object<-new("pepsum")
	object@filetype="inspect"
	#object@fraction<-input
	######read files, append the file name to each row######
	object@data<-read.delim(input[1],header=TRUE)
	source.col<-length(object@data)+1
	object@data[,source.col]<-input[i]
	if(length(input)>=2){
		for(i in 2:length(object@fraction)){
			u<-read.delim(input[i],header=TRUE)
			source.col<-length(u)+1
			u[,source.col]<-input[i]
			object@data<-merge(u,object@data,all=TRUE,sort=FALSE)
		}
	}
	print(source.col)
	if(index=="TRUE"){object<-indexfiles(object)}
	object@fraction<-levels(object@data[,1])
	return(object)
	}










#########################################################################
#######function to read protein prospecter text output from file	#######
#######	090811 	added the ability to run Pep.suM after file read
#######			currently sets object@fraction to 
#######			levels(object@data$Fraction)
############		OR make separate function to combine Pepsum objects
#########################################################################

read.prospector=function(dir=getwd(),input=list.files(),index="TRUE",filter="FALSE",summary="FALSE"){
	object<-new("pepsum")
	object@filetype="prospector"
	setwd(dir)
	input.l<-length(input)
	#for(i in 1:input.l)
	con<-file(description = input, blocking = FALSE, open = "r")
	lines<-readLines(con)
	###find the line before the table and add one to it
	####change this to add 2 if actually comparing searches
	skiplines<-grep("Search Name:",lines)+1
	skiplines<-max(skiplines)
	object@data<-read.delim(input,fill=TRUE,skip=skiplines,header=TRUE,row.names=NULL)
	close(con)
	object@fraction<-levels(object@data$Fraction)
	object@fraction<-object@fraction[object@fraction!=""]
	if(filter=="TRUE"){
		object<-filterfiles(object)
		}
	if(index=="TRUE"){object<-indexfiles(object)}
	if(summary=="TRUE"){
		object<-Pep.suM(object)
		print(object@summary)
		}
	return(object)
	}



#########################################################################
################          read pilot				#############
################		UPDATED:050411   ################################
#########################################################################

read.pilot=function(dir=getwd(), input=list.files(),conf=99,index="TRUE"){
	require(gdata)
	setwd(dir)
	####make new object and start filling it################
	object<-new("pepsum")
	object@filetype="pilot"
	object@filenames<-input
	######read files, append the file name to each row######
	object@data<-read.delim(input[1],header=TRUE)
	source.col<-length(object@data)+1
	object@data[,source.col]<-input[i]
	if(length(input)>=2){
		for(i in 2:length(object@fraction)){
			u<-read.delim(input[i],header=TRUE)
			source.col<-length(u)+1
			u[,source.col]<-input[i]
			object@data<-merge(u,object@data,all=TRUE,sort=FALSE)
		}
	}
	print(source.col)
	if(index=="TRUE"){object<-indexfiles(object)}
	return(object)
	}


#########################################################################
#############function to generate index of values above conf	#######
#########	runs very fast
#########################################################################
###return a matrix with column indexes that pass that "by"
##########
#####		12/21/2011 this works to index proteins by activation, accession, fractions, and pval
###############################################################

####  05192012	 fixed to add flexibility in finding the fdr column

indexfiles=function(object,by=list(pval=0.01,fraction=c(),proteins=c("P02769","P01966","P02070","P02662","P02663","P02666","Tryp_Pig","P00698","Q6B411"),activation=c("ETD","CID"))){
	
	if(object@filetype=="MSGFDB"){
		activations.vector<-c()
		if(length(by$activation)>=1){
			activations.vector<-object@data$FragMethod
			activation.levels<-levels(activations.vector)
			act.levels.length<-length(activation.levels)
			activations.length<-length(activations.vector)
			activation.index<-matrix(nrow=activations.length,ncol=act.levels.length,dimnames=list(c(),c(activation.levels)))
			for(i in 1:act.levels.length){
				indexes<-c(rep(0,activations.length))
				for(j in 1:activations.length){
					if(activations.vector[j]==activation.levels[i]){
						indexes[j]<-1
						}
					}
				activation.index[,i]<-indexes
				}
			print("indexed by activation")
			}
			### remove the option to index only certain fractions
		if(length(by$fraction)<1){
			fractions<-object@fraction
			frac.l<-length(fractions)
			data.l<-length(object@data[,"Peptide"])
			fracindex<-matrix(nrow=data.l,ncol=frac.l,dimnames=list(c(),c(fractions)))
			fractioncolumn<-length(object@data[1,])
			for(j in 1:frac.l){
				indexes<-c(rep(0,data.l))
				for(i in 1:data.l){
					if(object@data[i,fractioncolumn]==fractions[j]){
						indexes[i]<-1
						}
					}
				fracindex[,j]<-indexes
				}
			print("indexed by all fractions")
			}

		#if(length(by$fraction)>=1){
			#print("so far so good")
			#fractions<-by$fraction
			#indexes<-c()
			#frac.l<-length(fractions)
			#}
		if(length(by$pval>=1)){
			###make an index of peps below 1%
			pval<-by$pval
			fractions<-c()
			fdrcolumn<-length(object@data[1,])-1
			pepfdrs<-object@data[,fdrcolumn]
			pepfdrs.l<-length(pepfdrs)
			index<-c(rep(0,data.l))
			k<-1
			binaryfdr<-c(rep(0,data.l))
			for(i in 1:pepfdrs.l){
				if(pepfdrs[i]<=pval){
					index[k]=i
					binaryfdr[i]<-1
					k=k+1
					}
				}
			print("indexed by pvalue")
			print(by$pval)
			}
		proteins<-by$proteins
		if(length(by$proteins>=1)){
			proteins<-by$proteins
			pro.l<-length(proteins)
			pro.vec<-object@data$Protein
			promatindex<-matrix(nrow=data.l,ncol=pro.l,dimnames=list(c(),c(proteins)))
			for(j in 1:pro.l){
				proindex<-c(rep(0,data.l))
				for(i in 1:data.l){
					if(length(grep(pattern=proteins[j],x=pro.vec[i],ignore.case=TRUE))>=1){
						proindex[i]<-1
						}
					}
				promatindex[,j]<-proindex
				}
			print("indexed by user protein accessions")
			}
		if(length(fracindex)>=1 && length(indexes)>=1 && length(proteins)>=1 && length(activations.vector)>=1){
			object@index<-list(pepfdrs=index,fractions=fracindex,proteins=promatindex,binaryfdr=binaryfdr,activation=activation.index)
			}
		if(length(fracindex)>=1 && length(indexes)>=1 && length(proteins)==0 && length(activations.vector)==0){
			object@index<-list(pepfdrs=index,fractions=fracindex,binaryfdr=binaryfdr)
			}
		if(length(fracindex)>=1 && length(indexes)>=1 && length(proteins)>=1 && length(activations.vector)==0){
			object@index<-list(pepfdrs=index,fractions=fracindex,proteins=promatindex,binaryfdr=binaryfdr)
			}
		}
	if(object@filetype=="pilot"){
		###make an index of peps at conf
		vec<-object@data$Conf
		vec[is.na(vec)]<-0
		m<-length(object@data$Conf)
		indexes<-c()
		k<-1
		for(i in 1:m){
			if(vec[i]>=99){
				indexes[k]=i
				k=k+1
				}
			}
		length.indexes<-length(indexes)#####how many values are in the new array
		###prove it worked
		#for(i in 1:length.indexes){
			#print(object@data[indexes[i],])
			#}
		object@index<-list(conf=indexes)
		return(object)
		}
	
	if(object@filetype=="prospector"){
	####make an index of peps in each fraction
		fraclist<-as.character(object@data$Fraction)
		fraclist.l<-length(fraclist)
		frac.levs<-as.character(object@fraction)
		frac.levs.l<-length(frac.levs)
		fracindex<-matrix(nrow=fraclist.l,ncol=frac.levs.l,dimnames=list(c(),frac.levs))
		for(j in 1:frac.levs.l){
			indexes<-c()
			k<-1
			for(i in 1:fraclist.l){
				if(fraclist[i]==frac.levs[j])indexes[k]<-i
				if(fraclist[i]!=frac.levs[j])indexes[k]<-0
				k=k+1
				}
			#####set the new indexes vector to the j-th slot in object@index
			fracindex[,j]=indexes
			####top of the loop to (j+1)
			}
		object@index<-list(fraction=fracindex,pepfdrs=seq(from=1,to=length(object@data[,1]),by=1))
		}
	return(object)
	}




#########################################################################
#############function to generate index of values above conf	#######
#########	runs very fast
#########	081711  added the ability to count FDRs instead of pepFDRs
#########	need to ask Sangtae about which to use....
#########################################################################
################v2

#########################################################################
######put values for each fraction into matrix 				#######
######column names= fraction or Sample ID					#######
######put NA in when the value is not part of that level		#######
#########################################################################

filterfiles=function(object=object){
	####	this is my older and slower method to get fraction specific sequences
	####	within multiple runs 
	####	makes a vector to loop through and a vector to compare against
	####  
	if(object@filetype=="prospector"){
		fractions<-c(as.character(object@data$Fraction))
		fractions<-fractions[fractions!=""]
		levs<-levels(object@data$Fraction)
		levs<-levs[levs!=""]
		object@sequence<-matrix(nrow=length(fractions),ncol=length(levs),dimnames=list(c(),c(levs)))
		peps<-object@data$DBPeptide
		for(x in levs){
			for(i in 1:length(fractions)){
				if(fractions[i]==x){object@sequence[i,x]<-c(as.character(peps[i]))}
				else if(fractions[i]!=x){object@sequence[i,x]<-NA}
				}
			}
		return(object)
		}
	else if(object@filetype=="pilot"){
		if(is.null(object@index)){object<-indexfiles(object)}
		fractions<-c(as.character(object@data$V27))
		fractions<-fractions[fractions!=""]
		levs<-object@fraction
		object@sequence<-matrix(nrow=length(fractions),ncol=length(levs),dimnames=list(c(),c(levs)))
		peps<-object@data$Sequence
		for(x in levs){
			for(i in 1:length(object@index)){
				if(fractions[object@index[i]]==x){object@sequence[object@index[i],x]<-c(as.character(peps[object@index[i]]))}
				else if(fractions[object@index[i]]!=x){object@sequence[object@index[i],x]<-NA}
				}
			}
		return(object)
		}
	}


#########################################################################
########	index sequences
#########################################################################

indexseq=function(object=object){
	####	I don't remember how this is different than the above script
	####
	####	makes a vector to loop through and a vector to compare against
	####

	if(object@filetype=="prospector"){
		fractions<-c(as.character(object@data$Fraction))
		fractions<-fractions[fractions!=""]
		levs<-levels(object@data$Fraction)
		levs<-levs[levs!=""]
		object@sequence<-matrix(nrow=length(fractions),ncol=length(levs),dimnames=list(c(),c(levs)))
		for(x in levs){
			for(i in 1:length(fractions)){
				if(fractions[i]==x){object@sequence[i,x]<-c(as.character(object@data$DBPeptide[i]))}
				else if(fractions[i]!=x){object@sequence[i,x]<-NA}
				}
			}
		return(object)
		}
	else if(object@filetype=="pilot"){
		if(is.null(object@index)){object<-indexfiles(object)}
		
		fractions<-c(as.character(object@data$V27))
		fractions<-fractions[fractions!=""]
		levs<-object@fraction
		object@sequence<-matrix(nrow=length(fractions),ncol=length(levs),dimnames=list(c(),c(levs)))
		seq<-c(as.character(object@data$Sequence))
		for(x in levs){
			for(i in 1:length(object@index)){
				if(fractions[object@index[i]]==x){object@sequence[object@index[i],x]<-seq[object@index[i]]}
				else if(fractions[object@index[i]]!=x){object@sequence[object@index[i],x]<-NA}
				}
			}
		return(object)
		}
	}


#########################################################################
####parse peptide output from database search and count amino locations##
####outputs textplot, heatmap, and tab delim .txt file of observed P8-P8'

strReverse=function(x){
	sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
}


#########################################################################
######	works with pilot and prospector summaries
#########################################################################
######	Function to generate heatmaps of observed specificy
######	useful when:  you are studying a new protease
#########################################################################
#####		usage: 	cutsites(input=youPepsumObject)
#########################################################################

cutsites=function(object,length=6,len=6,unique=T,exclude=c(),cor="TRUE",hist=c("FALSE","P1")){
	require(gplots)
	title<-object@filetype
	AA.charlist<-c("A","S","T","G","V","C","N","L","I","M","P","Y","W","Q","F","D","E","H","K","R")
	######exclude unwanted values
	fractions<-object@fraction
	fractions=fractions[fractions!=""]
	####require(pepsum)
	###make an empty summary matrix################################
	filetot<-matrix(c(rep(0,240)),nrow=20,dimnames=list(AA.charlist,
	c("P6","P5","P4","P3","P2","P1","P1'","P2'","P3'","P4'","P5'","P6'")),
	ncol=12)
	#### section for MSGFDB type object####
	if(object@filetype=="MSGFDB"){
		### this function cleans peptides and returns significant seqs
		peps<-peptide.string.parse(object)
		### only peptides above/= length
		peps<-peps[nchar(peps)>=len]
		#hist(nchar(unique(peps)),breaks=c(6:40))
		if(unique==T){
			peps<-unique(peps)
			peps<-na.omit(peps)
			}
		rev.vec<-strReverse(peps)
		positions<-list()
		revpos<-list()
		for(i in 1:length){
			positions[[i]]<-substr(peps,i,i)
			revpos[[i]]<-substr(rev.vec,i,i)
			}
		
		### use the line below to get a multiple sequence alignment of peptides around
		### the site of cleavage
		msa<-rep(0, times=length(peps))
		for(i in 1:length(peps)){
			endpos<-nchar(peps[i])
			startpos<-endpos-len+1
			msa[i]<-paste(substr(peps[i],startpos,endpos),substr(peps[i],1,len),sep="")
			}

as.factor(msa)
?write.table
write.table(as.factor(msa), sep="\n",file="aLPMSA.txt",row.names=F,col.names=F,quote=F)

msa



		listed<-list(revpos[[6]],revpos[[5]],revpos[[4]],revpos[[3]],revpos[[2]],
			revpos[[1]],positions[[1]],positions[[2]],positions[[3]],positions[[4]],
			positions[[5]],positions[[6]])
		}

	if(object@filetype=="inspect"){
		seq.all<-object@data$Annotation
		seq.all<-unique(seq.all)
		seq.all<-as.character(seq.all)
		sized<-c()
		vec.l<-length(seq.all)
		for(i in 1:vec.l)
			if(nchar(seq.all[i])>=length+4)
				sized[i]<-seq.all[i]

		sized<-unique(sized)
		sized<-na.omit(sized)
		sized<-as.character(sized)
		rev.vec<-strReverse(sized)
		p1p<-substr(sized,3,3)
		p2p<-substr(sized,4,4)
		p3p<-substr(sized,5,5)
		p4p<-substr(sized,6,6)
		p5p<-substr(sized,7,7)
		p6p<-substr(sized,8,8)
		p7p<-substr(sized,9,9)
		p8p<-substr(sized,10,10)
		########reverse strings to find c-term P8-P1##########
		p1<-substr(rev.vec,3,3)
		p2<-substr(rev.vec,4,4)
		p3<-substr(rev.vec,5,5)
		p4<-substr(rev.vec,6,6)
		p5<-substr(rev.vec,7,7)
		p6<-substr(rev.vec,8,8)
		p7<-substr(rev.vec,9,9)
		p8<-substr(rev.vec,10,10)
		listed<-list(p8,p7,p6,p5,p4,p3,p2,p1,p1p,p2p,p3p,p4p,p5p,p6p,p7p,p8p)
		}
	if(object@filetype=="prospector"){
		#if(length(exclude)>=1){
		#######this works but doesnt change the output
		#######need to correlate it with the values of DBpep.all
		row.index<-c()
		for(i in 1:length(exclude))
			exclude[i]==object@fraction[j]	

		####	prev amino acids
		####	use these or not to increase the number of observations in cutsites
		#prevAA<-object@data$Prev.AA
		#nextAA<-object@data$Next.AA
		####	get list from object@data
		DBpep.all<-as.character(object@data$DB.Peptide)
		#DBpep.excluded<-DBpep.all[object@index[,1]]
		DBpep.all<-unique(DBpep.all)
		DBpep.all<-as.character(DBpep.all)						
		sized<-c()
		vec.l<-length(DBpep.all)
		for(i in 1:vec.l)
			if(nchar(DBpep.all[i])>=length) 
				sized[i]<-DBpep.all[i]
		########remove duplicates and NAs##########
		sized<-unique(sized)
		sized<-na.omit(sized)
		sized<-as.character(sized)
		rev.vec<-strReverse(sized)
		############REWRITE THIS as a LOOP over for(i in 1:length(listed[i]))		
		#####get vectors of single residues at P1'-P8'#####
		p1p<-substr(sized,1,1)
		p2p<-substr(sized,2,2)
		p3p<-substr(sized,3,3)
		p4p<-substr(sized,4,4)
		p5p<-substr(sized,5,5)
		p6p<-substr(sized,6,6)
		p7p<-substr(sized,7,7)
		p8p<-substr(sized,8,8)
		########reverse strings to find c-term P8-P1##########
		p1<-substr(rev.vec,1,1)
		p2<-substr(rev.vec,2,2)
		p3<-substr(rev.vec,3,3)
		p4<-substr(rev.vec,4,4)
		p5<-substr(rev.vec,5,5)
		p6<-substr(rev.vec,6,6)
		p7<-substr(rev.vec,7,7)
		p8<-substr(rev.vec,8,8)
		listed<-list(p8,p7,p6,p5,p4,p3,p2,p1,p1p,p2p,p3p,p4p,p5p,p6p,p7p,p8p)
		}
	if(object@filetype=="pilot"){
		####Look for fractions larger than 8 and 
		object@sequence<-c()
		sequences<-as.character(object@data$Sequence)
		sized<-c()
		k=1
		for(x in object@index){
		if(nchar(sequences[x])>8) 
			sized<-c(sized,x)
			}
		sized.l<-length(sized)
		for(i in 1:sized.l)sized[i]<-sequences[sized[i]]

		########remove duplicates and NAs##########
		sized<-unique(sized)
		sized<-na.omit(sized)
		sized<-as.character(sized)
		rev.vec<-strReverse(sized)
		############REWRITE THIS as a LOOP over for(i in 1:length(listed[i]))		
		#####get vectors of single residues at P1'-P8'#####
		p1p<-substr(sized,1,1)
		p2p<-substr(sized,2,2)
		p3p<-substr(sized,3,3)
		p4p<-substr(sized,4,4)
		p5p<-substr(sized,5,5)
		p6p<-substr(sized,6,6)
		p7p<-substr(sized,7,7)
		p8p<-substr(sized,8,8)
		########reverse strings to find c-term P8-P1##########
		p1<-substr(rev.vec,1,1)
		p2<-substr(rev.vec,2,2)
		p3<-substr(rev.vec,3,3)
		p4<-substr(rev.vec,4,4)
		p5<-substr(rev.vec,5,5)
		p6<-substr(rev.vec,6,6)
		p7<-substr(rev.vec,7,7)
		p8<-substr(rev.vec,8,8)
		listed<-list(p8,p7,p6,p5,p4,p3,p2,p1,p1p,p2p,p3p,p4p,p5p,p6p,p7p,p8p)
		}
		##########Make a matrix of P8-P8'###########################
	count<-0
	summary<-matrix(c(rep(0,240)),nrow=20,dimnames=list(
	c("A","S","T","G","V","C","N","L","I","M","P","Y","W","Q","F","D","E","H","K","R"),
	c("P6","P5","P4","P3","P2","P1","P1'","P2'","P3'","P4'","P5'","P6'")),
	ncol=12)
	for(x in listed){
      	count<-count+1
		summary[,count]<-c(
      	length(x[x=="A"]),length(x[x=="S"]),length(x[x=="T"]),
      	length(x[x=="G"]),length(x[x=="V"]),length(x[x=="C"]),
      	length(x[x=="N"]),length(x[x=="L"]),length(x[x=="I"]),
      	length(x[x=="M"]),length(x[x=="P"]),length(x[x=="Y"]),
      	length(x[x=="W"]),length(x[x=="Q"]),length(x[x=="F"]),
      	length(x[x=="D"]),length(x[x=="E"]),length(x[x=="H"]),
      	length(x[x=="K"]),length(x[x=="R"]))
		}
	filetot<-summary+filetot
	print(summary)
	print(filetot)
	object@totalheat<-filetot

	print(length(sized))
	##########This part corrects for observed column totals###
	if(cor=="TRUE"){
		dims<-dimnames(filetot)
		a.acids<-unlist(dims[1])
		cor.mat<-matrix(nrow=20,ncol=1)
		for(i in 1:length(cor.mat)){
			cor.mat[i]<-sum(filetot[i,])
			}
		for(j in 1:length(cor.mat)){
			for(i in 1:length(filetot[j,])){
				filetot[j,i]<-filetot[j,i]/cor.mat[j]
				}
			}
		
		rounded<-round(filetot,4)
		pct<-rounded*10000
		############HeatMap for corrected###################################
		jpeg(file=paste(title,"_cor",".jpeg",sep=""),quality=100,width = 480, height = 480, units = "px", pointsize = 12, 
		bg = "white")
		heatmap.2(filetot,main="Corrected for Amino Acid abundance",cellnote=pct,Rowv=NA,Colv=NA,trace=c("none"),margins=c(5,5), colsep=c(8),sepcol="blue",tracecol="blue",dendrogram = c("none"),ylab="Amino Acid",cexRow=2,cexCol=2,xlab = "site relative to cleavage",key=TRUE,keysize=4,density.info=c("density"),notecex=2.0,notecol="black",na.color=par("bg"))
		dev.off()
		############ratiotextplot
		#return(object)
		}

	if(cor=="TRUE"){
		dims<-dimnames(filetot)
		a.acids<-unlist(dims[1])
		cor.mat<-matrix(nrow=20,ncol=1)
		for(i in 1:length(cor.mat)){
			cor.mat[i]<-sum(filetot[i,])
			}
		for(j in 1:length(cor.mat)){
			for(i in 1:length(filetot[j,])){
				filetot[j,i]<-filetot[j,i]/cor.mat[j]
				}
			}
		rounded<-round(filetot,2)
		pct<-rounded*100
		############HeatMap for corrected###################################
		jpeg(file=paste(title,"_cor",".jpeg",sep=""),quality=100,width = 480, height = 480, units = "px", pointsize = 12, 
		bg = "white")
		heatmap.2(filetot,main="Corrected for Amino Acid abundance",cellnote=pct,Rowv=NA,Colv=NA,trace=c("none"),margins=c(5,5), colsep=c(8),sepcol="blue",tracecol="blue",dendrogram = c("none"),ylab="Amino Acid",cexRow=2,cexCol=2,xlab = "site relative to cleavage",key=TRUE,keysize=4,density.info=c("density"),notecex=2.0,notecol="black",na.color=par("bg"))
		dev.off()
		############ratiotextplot
		#return(object)
		}
	else if(cor=="FALSE"){
		############HeatMap for raw counts###################################
		jpeg(file=paste(title,"_raw",".jpeg",sep=""),width=900,height=900,units = "px", pointsize = 10, 
		quality = 100, bg = "white")
		heatmap.2(filetot,main="Raw counts: Observed Amino Acids",Rowv=NA,Colv=NA,trace=c("none"),colsep=c(8),sepcol="blue",ylab="Amino Acid",tracecol="blue",cexRow=2,cexCol=2,margins=c(5,5),xlab="site relative to cleavage",key=TRUE,keysize=4,density.info=c("density"),cellnote=filetot,notecex=1.5,notecol="black",na.color=par("bg"))
		dev.off()
		}
	############Total Summary###############################
	rounded<-round(filetot,3)
	write.table(rounded,file=paste("specificity matrix",".txt",sep=""),sep="\t")
	object@totalheat=filetot
	return(object)
	}







##################################################################################
##################################  082011
########			this function works with MSGFDB type pepsum objects
########	090411	adding prospecter to the filetypes this function can proc
########			added exclude=c() to remove fractions you don't want
########			works! needed a "" in the exclude default LOL
#####################################################################

charge.count=function(object,exclude=c("")){
	if(object@filetype=="MSGFDB"){
		###use object@index as the value of i during the fraction loop
		charges.all<-c()
		frac<-object@index$pepfdrs
		frac.l<-length(object@index$pepfdrs)
		listed<-c()
		for(x in object@fraction){
			peps<-c()
			charges<-c()
			lengths<-c()
			for(y in frac){
				if(object@data$V16[y]==x){
					#peps<-c(peps,as.character(object@data$Peptide[y]))
					charges<-c(charges,object@data$Charge[y])
					#lengths<-c(lengths,length(as.character(object@data$Peptide[y])))
					}
				}
			charges.all<-c(charges, charges.all)
			}
		charges.crunched<-c()
		for(i in 2:7){charges.crunched<-c(charges.crunched,length(charges.all[charges.all==i]))}
		for(i in 1:6){print(charges.crunched[i]/sum(charges.crunched))}
		pie(charges.crunched,labels=c(2,3,4,5,6,7),main="all charges observed")
		}
	if(object@filetype=="prospector"){
		charges.all<-c()
		peps.all<-c()
		lengths.all<-c()
		basic.residue.count<-c()
		pep.index<-c()
		fracs<-object@fraction[object@fraction!=exclude]
		##### make a vector index of peptides only in the vector fracs
		for(x in fracs){
			pep.index<-c(pep.index,object@index[,x])
			pep.index<-pep.index[!pep.index==0]
			}
		peps<-c()
		charges<-c()
		lengths<-c()
		for(y in pep.index){
			peps<-c(peps,as.character(object@data$DB.Peptide[y]))
			charges<-c(charges,object@data$z[y])
			}
		###section to count the bases in the peps
		lengths<-nchar(peps)
		peps<-strsplit(peps,split="")
		for(x in peps){
			basic.residue.count<-c(basic.residue.count,sum(length(x[x=="H"]),length(x[x=="K"]),length(x[x=="R"])))
			}
		length(basic.residue.count)
		plot(charges,basic.residue.count,main="correlation btw charge and basic residues")
		lm<-lm(charges~basic.residue.count)
		lm1<-lm(charges~basic.residue.count)
		anova(lm1)
		print(summary(lm))
		abline(b=coef(lm)[2],a=coef(lm)[1])
		
		}
	}




############################################################################
####		09/04/2011 working pep.sum script!!
#######		use with MSFGDB and prospector type objects
############################################################################	
##### 	09/08/11 	updated Pep.suM to return the peps as object@peps
######			and the summary as object@summary...
######			...Still need to fix the part about MSGFBD
############################################################################
######


Pep.suM=function(object){
	if(object@filetype=="prospector"){
		frac.l<-length(object@data$Fraction)
		fraclists<-as.character(object@data$Fraction)
		listed<-c()
		frac.levs<-levels(object@data$Fraction)
		for(x in frac.levs){
			peps<-c()
			charges<-c()
			lengths<-c()
			for(i in 1:frac.l){
				if(object@data$Fraction[i]==x){
					peps<-c(peps,as.character(object@data$DB.Peptide[i]))
					charges<-c(charges,object@data$z[i])
					lengths<-c(lengths,object@data$Length[i])
					}
				}
			listed<-c(listed,length(peps),length(unique(peps)),mean(lengths),mean(charges),sd(charges))
			}
		datmat<-matrix(listed,ncol=length(frac.levs),dimnames=list(c("# peptides","# unique peptides","ave. peptide length","ave. peptide charge","st. dev. charge"),c(frac.levs)))
		write.table(datmat,file="datmat.txt",sep="\t")
		object@summary<-datmat
		object@sequence<-peps
		return(object)
		}
	if(object@filetype=="MSGFDB"){
		frac.l<-length(object@fraction)
		fraclists<-as.character(object@fraction)
		listed<-c()
		for(x in fraclists){
			peptides<-c()
			charges<-c()
			lengths<-c()
			for(i in 1:frac.l){
				if(object@data$V16[i]==x){
					peptides<-c(peptides,as.character(object@data$Peptide[i]))
					charges<-c(charges,object@data$Charge[i])
					lengths<-c(lengths,length(object@data$Peptide[i]))
				}
			}
		}
	listed<-c(listed,x,length(peps),length(unique(peps)),mean(lengths),mean(charges),sd(charges))
	print(x)
	print(mean(charges))
	print(sd(charges))
	print(length(peps))
	print(length(unique(peps)))
	print(mean(lengths))
	
	datmat<-matrix(listed,nrow=6)
	write.table(datmat,file="datmat.txt",sep="\t")
	}
}


#########################################################################
######	Pep.meR 090811 kinda works....
#########################################################################
######	Function to build vector of peptides from both samples
######	useful when: comparing search output .txt files
#########################################################################
#####		usage:  peptide comparison
#########################################################################


Pep.meR=function(objects=c()){
	objects.l<-length(objects)
	peplist<-list()
	for(i in 1:objects.l){
		object<-objects[[i]]
		peplist[i]<-list(object@sequence)
		names
		}
	print(peplist)
	print(object@filetype)
	return(peplist)
	}

#########################################################################
######	Eic.calC
#########################################################################
######	Function to build vector of peptides from both samples
######	useful when: 
#########################################################################
#####		usage:  peptide comparison
#########################################################################

Eic.calC=function(object=pepsum){
	require(xcms)





	}


#########################################################################
######	below function doesnt work yet!!!!! 08/21/11
#########################################################################
######	Function to use cutsites frequency output and generate peptides
######	
#########################################################################
#####		usage:  to generate peptides by modelling lengths due given
######			observed frequency tables from cutsites()
#########################################################################


model.cut.proteome=function(object=Spombe, freq=object@totalheat[,8]){
	###setup objects
	freq.l<-length(freq)
	###only take the top 10 frequencies 
	freq.sorted<-sort(freq)
	###loop over each amino acid
	for(i in 1:20)
	
		rand.length<-object@AAratios[i,1]
		randint(
		#### loop over residues in freq
		for(x in freq){
			print(residues[1:n])
			peps.l<-length(peps[peps!=0])
			m<-1
			for(i in 1:peps.l){
        			pep<-unlist(strsplit(peps[i],split=""))
				pep.l<-length(pep)
				###loop over AAs in pro[i]
				k<-1
				for(j in 1:pep.l){
					if(pep[j]==residues[n]){
					##if the jth residue=residues
					###write pep 1:j
						peps2[m]<-c(paste(pep[k:j],collapse=""))
						k<-j+1
						m<-m+1
						}
					###if the jth residue == the last residue
					else if(j==pep.l){
						peps2[m]<-c(paste(pep[k:j],collapse=""))
						k<-j+1
						m<-m+1
						}
					}
				}
			peps<-peps2
			}
		peps<-peps[peps!=0]
		object@cleaved<-list(name<-peps)
		return(object)
	
		}



##########################################################################
######	This may or may not work...!
######	09/04/2011 works in parts
#########################################################################
######	Function to generate matrix where columns are the fixed residue
######	and the row = count of p2 observed
######  
#########################################################################


p2hist=function(object,cor=="TRUE"){
			DBpep.all<-unique(DBpep.all)
			DBpep.all<-as.character(DBpep.all)
			rev.vec.all<-strReverse(DBpep.all)
			p1.all<-substr(rev.vec.all,1,1)
			p2.all<-substr(rev.vec.all,2,2)
			p3.all<-substr(rev.vec.all,3,3)
			p1.all.1<-length(p1.all)
		p2.counts<-matrix(c(rep(0,400)),nrow=20,ncol=20,dimnames=list(AA.charlist,AA.charlist))
			
			for(j in 1:20){
			p2givenp1<-c()
			print(AA.charlist[j])
				for(i in 1:p1.all.1){
					if(p1.all[i]==AA.charlist[j]){p2givenp1<-c(p2givenp1,p2.all[i])}
					for(k in 1:20){
						p2.counts[k,j]<-c(length(p2givenp1[p2givenp1==AA.charlist[k]]))}
				}
			
			}
	if(cor=="TRUE"){
		cor.mat<-matrix(nrow=1,ncol=20)
		for(i in 1:length(cor.mat)){
			cor.mat[i]<-sum(p2.counts[,i])
			}
		for(j in 1:length(cor.mat)){
			for(i in 1:length(p2.counts[,j])){
				if(cor.mat[j]!=0)p2.counts[i,j]<-p2.counts[i,j]/cor.mat[j]
				}
			}
###write summary and .jpeg for for each file####
		#jpeg(file=paste(x,"plot.jpeg",sep=""),width=500,height=500,units="px",pointsize=10,
		#quality=100,bg="white")
		#barplot(summary[,4],names.arg=c(
		#"A","S","T","G","V","C","N","L","I","M","P","Y","W","Q","F","D","E","H","K","R"),
		#main=x,xlab="amino acid",ylab="# of observed peptides")
		#dev.off()
		#write.table(summary,file=paste("cut",conf,x,sep="_"),sep="\t")
			}


p3givenp1=function(object,cor="TRUE"){
	p3.counts<-matrix(c(rep(0,400)),nrow=20,ncol=20,dimnames=list(AA.charlist,AA.charlist))
	for(j in 1:20){
	p3givenp1<-c()
	print(AA.charlist[j])
		for(i in 1:p1.all.1){
			if(p1.all[i]==AA.charlist[j]){p3givenp1<-c(p3givenp1,p3.all[i])}
			for(k in 1:20){
				p3.counts[k,j]<-c(length(p3givenp1[p3givenp1==AA.charlist[k]]))}
		}
			
			}
###hist(p3.counts[,]
	if(cor=="TRUE"){
		cor.mat<-matrix(nrow=1,ncol=20)
		for(i in 1:length(cor.mat)){
			cor.mat[i]<-sum(p3.counts[,i])
			}
		for(j in 1:length(cor.mat)){
			for(i in 1:length(p3.counts[,j])){
				if(cor.mat[j]!=0)p3.counts[i,j]<-p3.counts[i,j]/cor.mat[j]
				}
			}
		}
	jpeg(file=paste("P3givenP1",".jpeg",sep=""),width=700,height=700,units = "px", pointsize = 10, 
		quality = 100, bg = "white")
		heatmap.2(p3.counts,main="Corrected for # of peptides with each P1",Rowv=NA,Colv=TRUE,trace=c("none"), dendrogram = c("column"),ylab="P3 frequency given X",cexRow=2,cexCol=2,xlab = "P1 residue",key=TRUE,keysize=3,notecex=2.0,notecol="black",na.color=par("bg"))
		dev.off()
	return(p3.counts)
	}

####make heatmap of P2 given P1

jpeg(file=paste("P2givenP1",".jpeg",sep=""),width=700,height=700,units = "px", pointsize = 10, 
		quality = 100, bg = "white")
		heatmap.2(p2.counts,main="Corrected for # of peptides with each P1",Rowv=NA,Colv=TRUE,trace=c("none"), dendrogram = c("column"),ylab="P2 frequency",cexRow=2,cexCol=2,xlab = "P1 residue",key=TRUE,keysize=3,notecex=2.0,notecol="black",na.color=par("bg"))
		dev.off()


}


#########################################################################
######	Pep.meR
#########################################################################
######	Function to build vector of peptides from both samples
######	useful when: 
#########################################################################
#####		usage:  peptide comparison
#########################################################################


Pep.meR=function(objects=list()){
	object.l<-length(objects)
	
	for(i in 1:length(objects))
		object=objects
	print(object1@filetype)
	print(object2@filetype)
	}


#########################################################################
######	below function doesnt work yet!!!!! 05/14/11
#########################################################################
######	Function to build distributable R package
######	useful when:  you want to publish
#########################################################################
#####		usage:  not sure.....
#########################################################################


#?package.skeleton
#package.skeleton(list=c("strReverse","object1","object2","filterfiles","indexfiles", "read.prospector","read.pilot"), name="pepsum",path="C:/")

#########################Make this into a package

#package.skeleton(name="Pepsum", path="C:/Users/Jesse/Documents/R/win-library/2.12",list=c("obj2","obj1"),code_files=c("C:/users/jesse/documents/pepsum/"),force=TRUE)


#########################################################################
######	Eic.calC
#########################################################################
######	Function to build vector of peptides from both samples
######	useful when: 
#########################################################################
#####		usage:  peptide comparison
#########################################################################

#################################################################################
################
#######################
######this function takes mass to charge and charge to compute monoisotopic peptide mass
################## usage:

getM=function(m.z,z){
	m = m.z*z
	M = m-z*1.007276466812
	return(M)
	}
#M<-getM(m.z=mz,z=z)


###this function permutes the monoisotopic mass to the one to five charge states
getMzVal=function(pepmass=M,chargemax=5){
		listcharge<-c()
	for(i in 1:chargemax){
		listcharge[i] = (pepmass+i*1.00727)/i
		}
	return(listcharge)
	}


#getMzVal(pepmass=1470.70)->mzvals

#####combine the above functions to make a matrix of values
mzExpander=function(pepsum=ltqetd,hist=FALSE){
	if(object@filetype=="prospector"){
		datmat<-object@data
		m.z<-datmat[,"m.z"]
		z <- datmat[,"z"]
		M<-getM(m.z,z)
		M.l<-length(M)
		predMzmat<-matrix(ncol=M.l,nrow=5,dimnames=list(c("+1","+2","+3","+4","+5"),c(M)))
		for(i in 1:M.l){
			mz.predicted<-getMzVal(pepmass=M[i])
			predMzmat[,i]<-mz.predicted
			}
		if(hist==TRUE){
			hist(M,breaks=40)
			}
		}
	if(object@filetype=="MSGFDB"){
		print("MSGFDB type pepsum objects not yet supported")
		}
	return(predMzmat)
	}
#expandedmzvalues<-mzExpander()





#####################################################
####   this gets the total spaces
##########
#####################################################################
########## usage: faTindex<-this(index<-c())


combinefractionindex=function(index=c()){
	m<-object@index$fraction
	row.length<-nrow(m)
	vecsum<-c(rep(0,row.length))
	for(x in index){
		vecsum<-vecsum+m[,x]
		}
	return(vecsum)
	}


###use this to add the columns of object@index$proteins


addmatrixcolumns=function(m=TEPms2indexed@index$proteins){
	row.length<-nrow(m)
	column.length<-ncol(m)
	vecsum<-c(rep(0,row.length))
	for(i in 1:column.length){
		vecsum<-vecsum+m[,i]
		}
	return(vecsum)
	}

##################################


###########################################################################
#####		change this such that it takes the fraction and protein names as input
#######################		old version
#############			input is a vector of index values 
#### determined by grep(object, pattern=length(criteria))
########################################
######  if pval is the only criteria, then use unique.peps instead



get.peps= function (input = object, index = pandpindexvals) {
	object = input
	if (object@filetype == "MSGFDB") {
		peps <- c()
		for (x in index) peps[x] <- as.character(object@data$Peptide[x])
		peps <- na.omit(peps)
		peps <- unique(peps)
		return(peps)
		}
	if (object@filetype != "MSGFDB") {
		print("wrong file type")
		}
	}


###########################################################################
#####		change this such that it takes the fraction and protein names as input
get.peps= function (input = object, proteins = c(), fractions = c()) {
	object = input
	### how many criteria are indexed?
	indexnames<-names(object@index)
	index.l<-length(indexnames)
	######
	if (object@filetype == "MSGFDB") {
	#########	first make a binary index of those criteria
		if(length(proteins)>=1){
			}
		
	#### grep for the length of the criteria

		peps <- c()
		for (x in index) peps[x] <- as.character(object@data$Peptide[x])
		peps <- na.omit(peps)
		peps <- unique(peps)
		return(peps)
		}
	if (object@filetype != "MSGFDB") {
		print("wrong file type")
		}
	}

####################################################################
###################################################################
###finally, filter gotten peps down to only those in all fractions
####################################################################
####################################################################
AA.masses<-list(A=71.03711,
	R=156.10111,
	N=114.04293,
	D=115.02694,
	C=103.00919,
	Q=128.05858,
	E=129.04259,
	G=57.02146,
	H=137.05891,
	I=113.08406,
	L=113.08406,
	K=128.09496,
	M=131.04049,
	F=147.06841,
	P=97.05276,
	S=87.03203,
	T=101.04768,
	W=186.07931,
	Y=163.06333,
	V=99.06841)
####################################################
################function to 
#######			calculate the monoisotopic mass from a charSeq
##########################################################################
###########################################################################
#####################################


calculatePepM=function(charSeq="TEDELQDKIHPF"){
	###table of exact amino acid masses
	AA.masses<-list(A=71.03711,
	R=156.10111,
	N=114.04293,
	D=115.02694,
	C=103.00919,
	Q=128.05858,
	E=129.04259,
	G=57.02146,
	H=137.05891,
	I=113.08406,
	L=113.08406,
	K=128.09496,
	M=131.04049,
	F=147.06841,
	P=97.05276,
	S=87.03203,
	T=101.04768,
	W=186.07931,
	Y=163.06333,
	V=99.06841)
	#######take the input character sequence	
	subsplt<-unlist(strsplit(charSeq,""))
	seq.len<-length(subsplt)
	sum=18.0105######water mass
	for(i in 1:seq.len){
		sum=sum+AA.masses[[subsplt[i]]]

		}
	print(sum)
	}

#calculatePepM()


#########################################################################
######	Does it work?   yes/no!!!!! [date]
#########################################################################
######	Function does: what????
######	useful when:  when?????
#########################################################################
#####		usage:  example usage???
#########################################################################

plotz=function(object=pepsum,type="m/z"){
	###start with just a m/z histogram#####
	mzvals<-c()
	scan.numz<-c()
	if(type=="m/z"){
		###setup vectors depending on the input filetype
		if(object@filetype=="MSGFDB"){
			for(i in 1:length(
			}
		###setup vectors depending on the input filetype
		if(object@filetype=="Prospector"){

			}

	
		
			}
		}
	###make a histogram of the values from above
	hist()
	}



#########################################################################
######	Does it work?   no!!!!! 
#########################################################################
######	Function does: writes a _syn.txt file for SpectrumLook
######	useful when:  you want to visualize your annotations 
#########################################################################
#####		usage:  write.syntxt()
#########################################################################



#syntxtnames<-colnames(syntxt)


write.syntxt=function(object=pepsum,mzXMLfiles=c("alp_sdc9.mzXML")){
	####first make a vector of the needed columns
	syntxtcolumns<-c("HitNum","ScanNum","ScanCount","ChargeState","MH","XCorr","DelCn","Sp","Reference","MultiProtein","Peptide","DelCn2","RankSp","RankXc","DelM","XcRatio","PassFilt","MScore","NumTrypticEnds")
	?data.frame
	if(object@filetype="MSGFDB"){
		fractionsindex<-c()
		for(x in object@data$pepfdrs)
			if(x==mzXMLfiles[1])
			object@index<-c(object@index,)
		###### make a data frame with the above names
		syntxt<-data.frame()
		colnames(syntxt)<-syntxtcolumns
			for(i in 1:length(object@index)){
				print(object@index)###check the contents of this slot
				}
			}


#########################################################################
######	Does it work?   no!!!!! 073111
#########################################################################
######	Function does: compiles peptides 
######	useful when:  compiling peptides
#########################################################################
#####		usage:  example usage???
#########################################################################



#########################################################################
######	Works for MSGFDB objects with object@pepvec 041212
#########################################################################
######	Function does: gets character-only strings from peptide vectors
######	useful when:  trying to compute % redundancy
#########################################################################
#####		usage:  parsedvector<-peptide.string.parse(object=read.MSGFDB())
#########################################################################

peptide.string.parse=function(object=alpm4){
	peptempvec<-object@pepvec
	peptempvec.l<-length(peptempvec)
	###loop over all the strings in pepvec
	peps<-c(rep(0,times=peptempvec.l))
	for(i in 1:peptempvec.l){
		### remove the first two and last two characters
		peps[i]<-substr(peptempvec[i],3,nchar(peptempvec[i])-2)
		###match the string starting with +/- followed by n integers, a period, and n more integers
		### and replace with nothing
		peps[i]<-gsub("(\\-|\\+)([0-9]+)(.)([0-9]+)", replacement="", peps[i]) 
		}
	return(peps)
	}

#peptide.string.parse()

#########################################################################
######	Works for any table of integers (e.g. from table()) 041312
#########################################################################
######	Function does: returns a normalized table 
######	useful when:  comparing raw distributions
#########################################################################
#####		usage:  newtable<-normalizetable(table=tryptic.lengths)
#########################################################################

normalizetable=function(table=tryptic.lengths){
	normalized<-c()
	for(i in 1:length(table))normalized<-c(normalized,table[i]/sum(table))
	return(normalized)
	}

#########################################################################
######	Works for any table of integers (e.g. from table()) 041312
#########################################################################
######	Function does: returns a normalized table 
######	useful when:  comparing raw distributions
#########################################################################
#####		usage:  newtable<-normalizetable(table=tryptic.lengths)
#########################################################################

chargecount=function(peptides=fapep){
	rkh<-c("R","K","H")
	basiccount<-c(rep(0,times=length(peptides)))
	for(i in 1:length(peptides)){
		string<-unlist(strsplit(fapep[i],split=""))
		basiccount[i]<-length(string[string==rkh[1]])+length(string[string==rkh[2]])+length(string[string==rkh[3]])
		
		}
	return(basiccount)
	}

#########################################################################
######	make the msgfdb peptides "nice"
#########################################################################
######	Function does: write a file with peptide sequences for PNNL
######	coverage summarizer program
######	useful when:  comparing raw distributions
#########################################################################
#####		usage:  newtable<-normalizetable(table=tryptic.lengths)
#########################################################################



alpboth
alp<-c(alpmat,alp2mat)
peptide.string.parse(aLP2)->alp2mat
peptide.string.parse(try2)->try2mat
peptide.string.parse(m190aboth)->m190amat

try<-c(try1mat,try2mat)
alptry<-c(alp,try)

write.table(as.factor(alptry),file="alptrybothmat.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names =FALSE)



