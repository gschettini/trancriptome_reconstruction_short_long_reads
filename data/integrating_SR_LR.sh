#!/bin/bash
source ./settings.sh

if [[ "$folder_b" != "" ]]; then
  echo "More than one set of samples were set or the default folders were set."
  source ./master/docker_script.sh
else
  echo "Only one set of samples was set."
  source ./master/docker_script_a.sh
fi
