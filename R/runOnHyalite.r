#Run chunks of code on hyalite from home computer
#Requires 'hyalite' to be set as an ssh alias


runHyalite <- function(jobId=NULL, objs=ls(envir=parent.frame()), packages=.packages(), updateMSUWC=T, user="jerad.hoy", cores=1, oneLine=F, overwrite=T){

	stp <- F
	code <- c()
	if(oneLine){
		print("Execute line to be sent to hyalite")
	} else {
		print("Execute lines to be sent to hyalite")
		print("Then type 'hyalite.send' to send code")
	}

	while(!stp){
		row <- readline()
		if(row == "hyalite.send"){
			stp <- T
		} else {
			code <- c(code, row)
		}
		if(oneLine == T){
			stp=T
		}
	}

	print("Processing objects")

	if(!exists(".jobs", envir=.GlobalEnv)){
		assign(".jobs", value=c(), envir=.GlobalEnv)
	}

	if(is.null(jobId)){
		jobId <- idmaker(1)
	}

	jobName <- paste0("r_", jobId)

	system(paste0("mkdir ", jobName))

	save(list=objs, file=paste0(jobName, "/rObjToSend.RData"), envir=parent.frame())
	writeLines(code, paste0(jobName, "/rScriptToRun.r"))

	writeLines(c("#!/bin/sh",
				 "#SBATCH -N 1",
				 paste0("#SBATCH -n ", cores),
				 paste0("#SBATCH -J ", jobId),
				 "#SBATCH --mail-user j24x165@montana.edu",
				 "#SBATCH --mail-type=ALL",
				 "#SBATCH --mail-type=END",
				 "#SBATCH  -o submitR.out",
				 "#SBATCH  -e submitR.err",
				 "#SBATCH -p priority",
				 "#SBATCH --no-requeue",
				 "",
				 paste0("srun -n 1 /home/", user, "/", jobName, "/", jobId, ".sh")),
			   paste0(jobName, "/sbatch_submitR.sh"))

	writeLines(c(
				if(updateMSUWC) paste0("devtools::install('/home/", user, "/msuwcRouting')"),
				paste0("fdir <- '/home/", user, "/", jobName, "/'"),
				paste0("library(", packages, ")"),
				"load(paste0(fdir, 'rObjToSend.RData'))",
				"runEnv <- new.env()",
				"source(paste0(fdir, 'rScriptToRun.r'), echo=T, local=runEnv, print.eval=T)",
				"save(list=ls(runEnv), file=paste0(fdir, 'rOutput.RData'), envir=runEnv)"),
			paste0(jobName, "/submitR.r"))

	writeLines(c(
				"#!/bin/sh",
				paste0("rm /home/", user, "/", jobName, "/job.out"),
				paste0("rm /home/", user, "/", jobName, "/rOutput.RData"),
				paste0("R --vanilla < /home/", user, "/", jobName, "/submitR.r &> /home/", user, "/", jobName, "/job.out")),
			paste0(jobName, "/", jobId, ".sh"))

	if(!(jobId %in% get(".jobs", envir=.GlobalEnv))){

		assign(".jobs", value=c(get(".jobs", envir=.GlobalEnv), jobId), envir=.GlobalEnv)
		print(paste0("Job ", jobId, " added to .jobs"))

	} else if(overwrite){
		print(paste0("Job ", jobId, " added to .jobs again!"))

	} else {
		response <- readline(paste0("Job '", jobId, "'already submitted, enter 'y' to overwrite?"))
		if(response != "y"){
			stop(paste0("Job '", jobId, "' submission cancelled"))
		}
	}

	if(updateMSUWC) system(paste0("rsync -avrP msuwcRouting hyalite:/home/", user, "/"))

	system(paste0("rsync -avrP ", jobName, " hyalite:/home/", user, "/"))

	system(paste0("rm -rf ", jobName))


	system(paste0("ssh hyalite 'chmod 555 /home/", user, "/", jobName, "/*'"))
	system(paste0("ssh hyalite 'sh /home/", user, "/", jobName, "/sbatch_submitR.sh' &"))

	print(paste0("Job ", jobId, " submitted to hyalite"))
	closeAllConnections()
}


getHyalite <- function(jList=.jobs, user="jerad.hoy", removeAfterLoad=F, showJobs=T){
	successes <- c()
	for(job in 1:length(jList)){

	    dat <- tryCatch(load(pipe(paste0(
					"ssh hyalite 'cat /home/", user, "/r_", jList[job], "/rOutput.RData'")),
					envir=.GlobalEnv, verbose=T), error=function(e){NULL})

		if(is.null(dat)){
			print(paste0("Job ", jList[job], " may not be finished or does not exist"))
		} else {
			print(paste0("Data from job ", jList[job], " loaded."))
			successes <- c(successes, job)
			if(removeAfterLoad){
				system(paste0("ssh hyalite 'rm -rf /home/", user, "/r_", jList[job], "'"))
			}
		}
	}

	print("Press a key to continue")
	invisible(readline())
	print(system("ssh hyalite sacct"))
	closeAllConnections()
	try(.jobs <<- .jobs[-c(successes)])
}


idmaker <- function(x){
	max.val = x*10000
	count <- nchar(as.character(max.val))

	# find out how many 'numbers' each ID will have after the letter
	size <- paste("%0",count,"d",sep="")
	# set the variable to be fed into 'sprintf' to ensure we have leading 0's
	lets <- toupper(sample(letters,x, replace=T))
	# randomising the letters
	nums <- sprintf(size,sample(1:max.val)[1])
	# randominsing the numbers, and ensuing they all have the same number of characters
	ids <- paste(lets,nums,sep="")
	# joining them together
	return(ids)
}
