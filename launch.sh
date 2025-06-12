#!/bin/bash

# For more info on how to set up this launch file, please refer to the README.md file in the github repository.

#====================================================================================#
# General setting, please set these variables to your own values                     #
# These are the settings that will be used by all scripts in the repository          #
#====================================================================================#

## Email notification settings
export email_notification="NO"                                  # either "YES" or "NO", in caps, If set to "YES", scripts in the repository will have line 5-6 changed to include your email address before submitting jobs.
export email_address="z5205618@ad.unsw.edu.au"                  # Set to your email address if you want notifications, otherwise leave it empty.
#export PROJECT="xl04"                                          # Should already be automatically set when you log in, but uncomment and set this manually if it's not automatically set.
## The path to the uniprot sprot and trembl made with diamond
export uniprot_sprot="/g/data/xl04/jc4878/database/uniprot_sprot.diamond.db.dmnd"        # Path to the uniprot_sprot diamond database, this should be the path to the diamond database file you created with the uniprot_sprot database. (leave unchanged if you have access to my directory)
export uniprot_trembl="/g/data/xl04/jc4878/database/uniprot_trembl.diamond.db.dmnd"      # Path to the uniprot_trembl diamond database, this should be the path to the diamond database file you created with the uniprot_trembl database. (leave unchanged if you have access to my directory)
## The path to your copy of the repository
export repository_path="/g/data/xl04/jc4878/github/Annotation_AusARG_Simplified"


#====================================================================================#
# Which pipeline are you running?     Only 1 can be launched at a time               #
#====================================================================================#

## P1 - runs script1, script2 then script3
export P1="NO"                      # either "YES" or "NO", in caps
## P2 - runs script1, then script3
export P2="NO"                      # either "YES" or "NO", in caps       

# Pipeline specific settings
## Set these if you are running P1
export P1_workingdir="/g/data/xl04/jc4878/workingdir"               # Do not put trailing "/" at the end, this is also the output directory
export P1_TRANSLATION_OUTPUT="/path/to/blastxtranslation/output"    # Do not put trailing "/" at the end, this should be set to the output directory from the blastxtranslation tool
export P1_GENOME="/path/to/genome.fasta.masked"                     # Path to your genome file, this should be the soft-masked genome file.
export P1_species="species_name"                                    # Species name, no space, can use prefixes like "TilRug" or full name "Tiliqua_rugosa"

## Set these if you are running P2
export P2_workingdir="/g/data/xl04/jc4878/workingdir"               # Do not put trailing "/" at the end, this is also the output directory
export P2_TRANSLATION_OUTPUT="/path/to/blastxtranslation/output"    # Do not put trailing "/" at the end, this should be set to the output directory from the blastxtranslation tool
export P2_GENOME="/path/to/genome.fasta.masked"                     # Path to your genome file, this should be the soft-masked genome file.
export P2_species="species_name"                                    # Species name, no space, can use prefixes like "TilRug" or full name "Tiliqua_rugosa"
export P2_AUGUSTUS_CONFIG="${P1_workingdir}/Augustus/config"        # Path to your augustus config folder, if you ran P1 previously for same species and so ${P1_workingdir} is also set above, then leave this as "${P1_workingdir}/Augustus/config", the script will automatically substitute the ${P1_workingdir} part to what you set above, you can also point this manually like "/path/to/Augustus/config"
export P2_DEPENDENCY_FOR_3="PBS_JOBID"                              # OPTIONAL, if this is set to anything other than "PBS_JOBID", then script3 in P2 will wait for the specified job to finish (in addition to script1) before starting. Leaving it as "PBS_JOBID" (default) will make script3 only wait for script1 to finish before starting.
## ^ P2_DEPENDENCY_FOR_3 is usually set to the job ID of script2 in P1 (if you are running P1 at the same time)


# Things to note
## script1 now runs minimap2 in parallel, 12 concurrent jobs at any time each using 1 thread, each job have peak RSS that varies from 10GB to 15GB, we request 190GB for the job which should be enough at any one time, but if you notice any memory problem in the PBS log, then simply decrease the ncpus to 6 or 8 and this should fix the problem.







#                              End of Configuration                              #
#                           Do not edit below this line                          #
#--------------------------------------------------------------------------------#

## Checking if anything is set to "YES" at all, if not, exit the script
if [ "$P1" = "YES" ] || [ "$P2" = "YES" ]; then
    echo -e "//===========================\\\\\\\\\n||          WELCOME          ||\n\\\\\\===========================//"
    export current_time=$(date)
    echo -e "${current_time}\n"
else
    echo -e "//===========================\\\\\\\\\n||          WELCOME          ||\n\\\\\\===========================//"
    export current_time=$(date)
    messages=(
    "[LOG] No tasks found. Taking a well-deserved nap... zzz"
    "[LOG] No tasks detected. System entering idle mode..."
    "[LOG] No tasks available. Initiating self-destruct sequence... Just kidding!"
    "[LOG] No tasks to process. Time for a coffee break!"
    "[LOG] No tasks detected. System is now in power-saving mode..."
    "[LOG] No tasks found. Enjoying a virtual vacation..."
    "[LOG] No tasks available. Playing solitaire until further notice..."
    "[LOG] No quests accepted. Returning to the main menu..."
    "[LOG] No tasks detected. Summoning cat videos instead..."
    )
    random_index=$((RANDOM % ${#messages[@]}))
    echo -e "${current_time}\n\n[LOG] Nothing set to 'YES'. Now exiting...\n${messages[$random_index]}"
    exit 0
fi

## Checking if only 1 pipeline is selected
selected_pipelines=0
[ "$P1" = "YES" ] && selected_pipelines=$((selected_pipelines+1))
[ "$P2" = "YES" ] && selected_pipelines=$((selected_pipelines+1))
echo -e "[LOG] Checking if you are running only one pipeline"
if [ "$selected_pipelines" -gt 1 ]; then
    echo "Error: Multiple pipelines set to 'YES'. Please select only one to launch at a time."
    exit 1
elif [ "$selected_pipelines" -eq 1 ]; then
    echo "[LOG] A pipeline detected, now preparing to launch..."
    exit 1
fi

## Adding email to scripts if email_notification is set to "YES"
if [ "$email_notification" = "YES" ]; then
    echo -e "[LOG] \$email_notification set to: YES, now adding email to scripts"
    sed -i "5s|.*|#PBS -M ${email_address}|" ${repository_path}/scripts/1GenerateTrainingGene.sh
    sed -i "5s|.*|#PBS -M ${email_address}|" ${repository_path}/scripts/2TrainingAugustus.sh
    sed -i "5s|.*|#PBS -M ${email_address}|" ${repository_path}/scripts/3RunningAugustus.sh
    sed -i "6s|.*|#PBS -m ae|" ${repository_path}/scripts/1GenerateTrainingGene.sh
    sed -i "6s|.*|#PBS -m ae|" ${repository_path}/scripts/2TrainingAugustus.sh
    sed -i "6s|.*|#PBS -m ae|" ${repository_path}/scripts/3RunningAugustus.sh
elif [ "$email_notification" = "NO" ]; then
    echo -e "[LOG] \$email_notification set to: NO, now removing email from scripts"
    sed -i "5s|.*||" ${repository_path}/scripts/1GenerateTrainingGene.sh
    sed -i "5s|.*||" ${repository_path}/scripts/2TrainingAugustus.sh
    sed -i "5s|.*||" ${repository_path}/scripts/3RunningAugustus.sh
    sed -i "6s|.*||" ${repository_path}/scripts/1GenerateTrainingGene.sh
    sed -i "6s|.*||" ${repository_path}/scripts/2TrainingAugustus.sh
    sed -i "6s|.*||" ${repository_path}/scripts/3RunningAugustus.sh
fi



#launching pipeline
if [ "$P1" = "YES" ]; then
    echo -e "[LOG] \$P1 set to: YES, now launching P1"
    mkdir -p ${P1_workingdir}
    mkdir -p ${P1_workingdir}/log
    mkdir -p ${P1_workingdir}/log/launch
    export P1_1_JOBID=$(qsub -P ${PROJECT} -o ${P1_workingdir}/log/PBS/P1_1GenerateTrainingGene.OU -v repository_path=${repository_path},workingdir=${P1_workingdir},TRANSLATION_OUTPUT=${P1_TRANSLATION_OUTPUT},genome=${P1_GENOME} ${repository_path}/scripts/1GenerateTrainingGene.sh)
    export P1_2_JOBID=$(qsub -P ${PROJECT} -W depend=afterok:${P1_1_JOBID} -o ${P1_workingdir}/log/PBS/P1_2TrainingAugustus.OU -v species=${P1_species},workingdir=${P1_workingdir},genome=${P1_GENOME} ${repository_path}/scripts/2TrainingAugustus.sh)
    qsub -P ${PROJECT} -W depend=afterok:${P1_2_JOBID} -o ${P1_workingdir}/log/PBS/P1_3RunningAugustus.OU -v species=${P1_species},workingdir=${P1_workingdir},genome=${P1_GENOME},AUGUSTUS_CONFIG=${P1_workingdir}/Augustus/config,CUSTOM_EXTRINSIC=${repository_path}/extrinsic.MPE_modded.cfg,uniprot_sprot=${uniprot_sprot},uniprot_trembl=${uniprot_trembl} ${repository_path}/scripts/3RunningAugustus.sh
    head -n43 "$0" > ${P1_workingdir}/log/launch/P1.launch.settings
    exit 0
elif [ "$P2" = "YES" ]; then
    echo -e "[LOG] \$P2 set to: YES, now launching P2"
    mkdir -p ${P2_workingdir}
    mkdir -p ${P2_workingdir}/log
    mkdir -p ${P2_workingdir}/log/launch
    export P2_1_JOBID=$(qsub -P ${PROJECT} -o ${P2_workingdir}/log/PBS/P2_1GenerateTrainingGene.OU -v repository_path=${repository_path},workingdir=${P2_workingdir},TRANSLATION_OUTPUT=${P2_TRANSLATION_OUTPUT},genome=${P2_GENOME} ${repository_path}/scripts/1GenerateTrainingGene.sh)
    if [ "$P2_DEPENDENCY_FOR_3" = "PBS_JOBID" ]; then
        qsub -P ${PROJECT} -W depend=afterok:${P2_1_JOBID} -o ${P2_workingdir}/log/PBS/P2_3RunningAugustus.OU -v species=${P2_species},workingdir=${P2_workingdir},genome=${P2_GENOME},AUGUSTUS_CONFIG=${P2_AUGUSTUS_CONFIG},CUSTOM_EXTRINSIC=${repository_path}/extrinsic.MPE_modded.cfg,uniprot_sprot=${uniprot_sprot},uniprot_trembl=${uniprot_trembl} ${repository_path}/scripts/3RunningAugustus.sh
    else
        qsub -P ${PROJECT} -W depend=afterok:${P2_1_JOBID}:${P2_DEPENDENCY_FOR_3} -o ${P2_workingdir}/log/PBS/P2_3RunningAugustus.OU -v species=${P2_species},workingdir=${P2_workingdir},genome=${P2_GENOME},AUGUSTUS_CONFIG=${P2_AUGUSTUS_CONFIG},CUSTOM_EXTRINSIC=${repository_path}/extrinsic.MPE_modded.cfg,uniprot_sprot=${uniprot_sprot},uniprot_trembl=${uniprot_trembl} ${repository_path}/scripts/3RunningAugustus.sh
    fi
    head -n43 "$0" > ${P2_workingdir}/log/launch/P2.launch.settings
    exit 0
fi

