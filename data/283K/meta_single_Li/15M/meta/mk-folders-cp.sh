#!/bin/bash
for ((i=1; i<=10; i++));
do
   dir_name="$(printf '%02d' $i)_IDNR"
   # dir_name="${i}_IDNR"
   mkdir $dir_name
   cp ./base_inputs/* $dir_name
done

