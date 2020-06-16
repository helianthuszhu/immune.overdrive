#####clean the clinical information
clin.data=data.frame(TCGA_subtypes=clin.CRC.full.sel$TCGA_subtypes,
                     mole_subtype=clin.CRC.full.sel$mole_subtype,
                     age=clin.CRC.full.sel$age_at_initial_pathologic_diagnosis,
                     gender=clin.CRC.full.sel$gender.y,
                     location=clin.CRC.full.sel$location,
                     hypermutation=clin.CRC.full.sel$hypermutation.x,
                     lymphatic_invasion=clin.CRC.full.sel$lymphatic_invasion.y,
                     pathologic_stage=clin.CRC.full.sel$pathologic_stage,
                     CMS.cluster=clin.CRC.full.sel$cms_label,
                     treatment_success=clin.CRC.full.sel$followup_treatment_success,
                     colon_polyps_present=clin.CRC.full.sel$colon_polyps_present,
                     histtory.NDJ.treatment=clin.CRC.full.sel$history_of_neoadjuvant_treatment,
                     person_neoplasm_cancer_status=clin.CRC.full.sel$person_neoplasm_cancer_status,
                     pathologic_M=clin.CRC.full.sel$pathologic_M,
                     pathologic_N=clin.CRC.full.sel$pathologic_N,
                     pathologic_T=clin.CRC.full.sel$pathologic_T,
                     cimp1=clin.CRC.full.sel$cimp,
                     cimpNature=clin.CRC.full.sel$CIMP,
                     MLH1_silencing=clin.CRC.full.sel$AWG_MLH1_silencing
                     )
rownames(clin.data)=rownames(clin.CRC.full.sel)
head(clin.data)
#stage
table(clin.data$pathologic_stage)
clin.data$pathologic_stage=gsub(" ", ".", clin.data$pathologic_stage)
clin.data$pathologic_stage=ifelse(grepl("Stage",clin.data$pathologic_stage),
                                  clin.data$pathologic_stage,NA)
library(DataCombine)
Replaces <- data.frame(from = names(table(clin.data$pathologic_stage)), 
                       to = c('stageI','stageI','stageII','stageII','stageII','stageII',
                              'stageIII','stageIII','stageIII','stageIII',
                              'stageIV','stageIV','stageIV'
                       ))
clin.data$AJCC.stage<- FindReplace(data = clin.data, Var = "pathologic_stage", replaceData = Replaces,
                                             from = "from", to = "to",vector=TRUE)
clin.data$AJCC.stage.bi=ifelse(clin.data$AJCC.stage=="stage1"|clin.data$AJCC.stage=="stage2","I&II","III&IV")
table(clin.data$AJCC.stage)
clin.data.CRC.selected=clin.data
save(clin.data.CRC.selected,file="clin.data.CRC.selected.RData")
#######
colnames(clin.CRC.full.sel)
