#!/bin/bash

# Some common functionality to avoid repeating it in each qsub

# Create logs folder in current working directory if not present
mkdir -p logs

# Set datafolder
datafolder="/projectnb/paxlab/EnhancerDiscovery/data"

# Used for default command-line arguments, e.g. `default $1 foo` will return the first argument, if it exists, otherwise foo.
function default() {
  if [[ -z "$1" ]];
    then echo $2
  else
    echo $1
  fi
}

# Start log entry with pre-run information. No arguments needed.
function prerun() {
  echo "=========================================================="
  echo "Starting on       : $(date)"
  echo "Run by            : $(whoami)"
  echo "Running on node   : $(hostname)"
  echo "Current directory : $(pwd)"
  echo "Current job ID    : $JOB_ID"
  echo "Current job name  : $JOB_NAME"
  echo "Task index number : $SGE_TASK_ID"
  echo "=========================================================="
}

# Use this function to log terminal commands, e.g. `run "echo foo"` will log the command being run, and then run the command.
# Helpful for contextualizing log outputs and errors in multi-step qsubs.
function run() {
  echo "Running command: '$1'"
  eval $1
}

# End log entry with post-run information. No arguments needed.
function postrun() {
  echo "=========================================================="
  echo "Finished on       : $(date)"
  echo "=========================================================="
}