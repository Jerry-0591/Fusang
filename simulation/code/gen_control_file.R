#!/usr/bin/env Rscript

### ARGUMENTS: N taxa N_sims Aln_length

args = commandArgs(trailingOnly=TRUE)
library('phangorn')
library('MCMCpack')
library('dplyr')
library('scales')
options(scipen=999)
#Model block generating function
model_gen=function(modelset,file,max_indel_length,indel_substitution_rate_lower_bound,indel_substitution_rate_upper_bound)
{
  len=max_indel_length
  model_id=1
  models_selected=c()
  for (model in modelset)
  {
    model_orig=model
    #Invariant sites Unif
    I=runif(1,0,1)
    A=runif(1,0,5)
    #Nucl proportions DIRICHLET 
    options(digits=5)
    Pi=format(rdirichlet(1, alpha=c(5,5,5,5)))
    #IndelRate = format(runif(1,indel_substitution_rate_lower_bound,indel_substitution_rate_upper_bound))

    models_selected=c(models_selected,paste(model,'Model',model_id,sep = ''))
    write(paste('\n[MODEL] ',model,'Model',model_id,sep = ''),file,append=T)
    options(digits=2)
    if (model %in% c('HKY','K80')){
      model=paste(c(model,' ',format(runif(1,0,3))),sep = '')
    } else if (model == 'TrN'){
      model=paste(c(model,' ',format(runif(2,0,3))),sep = '')
    } else if (model %in% c('TIM' ,'TIMef')){
      model=paste(c(model,' ',format(runif(3,0,3))),sep = '')
    } else if (model == 'TVM'){
      model=paste(c(model,' ',format(runif(4,0,3))),sep = '')
    } else if (model %in% c('SYM','GTR')){
      model=paste(c(model,' ',format(runif(5,0,3))),sep = '')
    } else if (model == 'UNREST'){
      model=paste(c(model,' ',format(runif(11,0,3))),sep = '')
    } else {
      model=model
    }
    model_id=model_id+1
    write(paste(' [submodel] ',paste(model,collapse=' '),'\n [rates] ',I,' ',A,' 0','\n [indelmodel] POW 1.5 ', paste(len,collapse=' '), '\n [indelrate] ', paste(runif(1,indel_substitution_rate_lower_bound,indel_substitution_rate_upper_bound),collapse='')),file,append=T)
    if (model_orig %in% c('F81','HKY','TrN','TIM','TVM','GTR'))
    {
      write(paste(' [statefreq]',paste(Pi,collapse=' ')),file,append=T)
    }
  }
  return(models_selected)
}

indelib_gen=function(n_taxa,n_sim,aln_length_upper_bound,aln_length_lower_bound,indel_substitution_rate_lower_bound,indel_substitution_rate_upper_bound,max_indel_length) # n_sim = number of simulations per topology
{
  
  write(paste('[TYPE] NUCLEOTIDE 2\n[SETTINGS]\n [output] FASTA\n [randomseed] ',round(runif(1,1,100000))),'control.txt')
  n_datasets=n_sim
  
  #Set MODEL block
  modelset=sample(c('JC','TIM','TIMef','GTR','UNREST'),n_datasets,replace=T)
  MODEL=model_gen(modelset,'control.txt',max_indel_length,indel_substitution_rate_lower_bound,indel_substitution_rate_upper_bound)
  
  #Set TREE block
  ID_TREE=paste(rep("t_sim",times=n_datasets),1:n_datasets,sep="")
  print("Newick")
  NEWICK=read.csv('../label_file/newick.csv',header=TRUE)
  NEWICK=NEWICK[,2]
  print("Done newick")
  write.table(data.frame('[TREE]',ID_TREE,NEWICK),'control.txt',append=T,quote=F,row.names=F,col.names =F)
  
  #Set PARTITIONS block
  PNAME=paste("p",1:n_datasets,sep="")
  write.table(data.frame('[PARTITIONS]',PNAME,"[",ID_TREE,MODEL,round(runif(n_sim,aln_length_lower_bound,aln_length_upper_bound)),"]"),'control.txt',append=T,quote=F,row.names=F,col.names =F)
  
  #Set EVOLVE block
  write('[EVOLVE]','control.txt',append=T)
  write.table(data.frame(PNAME,1,apply(data.frame(ID_TREE),1,paste,collapse="")),'control.txt',append=T,quote=F,row.names=F,col.names =F)
}

indelib_gen(as.numeric(args[1]),as.numeric(args[2]),as.numeric(args[3]), as.numeric(args[4]), as.numeric(args[5]), as.numeric(args[6]), as.numeric(args[7]))
