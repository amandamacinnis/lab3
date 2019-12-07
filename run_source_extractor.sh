#! /bin/bash -u

for file in $(ls -1 data/*.new)
do
   sex ${file} -c source_extractor_config_files/${file:5:-4}.sex -CATALOG_NAME source_extractor_output/${file:5:-4}.cat
done
