##########deal with KIRC clinal data
######
clin.KIRC=read.table("8.KIRC/KIRC_clinicalMatrix",header = T,sep = "\t")
head(clin.KIRC)
clin.KIRC.sel=data.frame(sampleID=clin.KIRC$sampleID,
                         age=clin.KIRC$age_at_initial_pathologic_diagnosis,
                         followup_treatment_success=clin.KIRC$followup_treatment_success,
                         gender=clin.KIRC$gender,
                         history_of_neoadjuvant_treatment=clin.KIRC$history_of_neoadjuvant_treatment,
                         laterality=clin.KIRC$laterality,
                         intermediate_dimension=clin.KIRC$intermediate_dimension,
                         karnofsky_performance_score=clin.KIRC$karnofsky_performance_score,
                         longest_dimension=clin.KIRC$longest_dimension,
                         neoplasm_histologic_grade=clin.KIRC$neoplasm_histologic_grade,
                         number_pack_years_smoked=clin.KIRC$number_pack_years_smoked,
                         pathologic_stage=clin.KIRC$pathologic_stage,
                         shortest_dimension=clin.KIRC$shortest_dimension,
                         id1=clin.KIRC$bcr_patient_barcode
                         )
rownames(clin.KIRC.sel)=clin.KIRC.sel$sampleID
clin.KIRC.sel$id2=substr(rownames(clin.KIRC.sel),14,15)
clin.KIRC.sel=subset(clin.KIRC.sel, id2=="01")
table(clin.KIRC.sel$id2)
head(clin.KIRC.sel)
#####molecular information
mole.KIRC=read.table("8.KIRC/clin/mole.subtype.snv.KIRC.txt",header = T,sep = "\t")
colnames(mole.KIRC)[1]="id1"
head(mole.KIRC)
#####
mole.KIRC.2=read.table("8.KIRC/clin/mole.cluster.levels.txt",header = T,sep = "\t")
mole.KIRC.2$id1=substr(mole.KIRC.2$Patient.tumor,1,12)
head(mole.KIRC.2)
dim(mole.KIRC.2)
length(unique(mole.KIRC.2$id1))
#combine two molecular infor
mole.KIRC.cb=merge(mole.KIRC, mole.KIRC.2, by="id1",all=TRUE)
head(mole.KIRC.cb)
table(mole.KIRC.cb$final.subtype.call.panRCC)
############
###########combine all
clin.KIRC.sel.full=merge(clin.KIRC.sel, mole.KIRC.cb, by="id1", all.x=TRUE)
rownames(clin.KIRC.sel.full)=clin.KIRC.sel.full$sampleID
head(clin.KIRC.sel.full)
dim(clin.KIRC.sel.full)
###########
save(clin.KIRC.sel.full,file="8.KIRC/clin/clin.KIRC.sel.full.537s.RData")
