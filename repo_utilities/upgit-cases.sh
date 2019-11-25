#!/bin/bash

# set test cases  directory
cases_dir="../cases"

for entry in "$cases_dir"/*
do
  echo "Updating gitignore in $entry..."
  cp gitignore.org $entry/.gitignore
done

echo "All done!"


