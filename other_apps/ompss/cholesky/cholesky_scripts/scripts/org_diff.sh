#!/bin/bash

# walk through all files in the current dir (and subdirs)
# and compare them with other files, showing percentage
# of duplication.

# which type files to compare?
# (wouldn't make sense to compare binary formats)
ext="txt"

# support filenames with spaces:
IFS=$(echo -en "\n\b")

working_dir="$PWD"
working_dir_name=$(echo $working_dir | sed 's|.*/||')
all_files="$working_dir/../$working_dir_name-filelist.txt"
remaining_files="$working_dir/../$working_dir_name-remaining.txt"

# get information about files:
find -type f -print0 | xargs -0 stat -c "%s %n" | grep -v "/\." | \
     grep "\.$ext" | sort -nr > $all_files

cp $all_files $remaining_files

while read string; do
    fileA="tasks1.txt"
#$(echo $string | sed 's/.[^.]*\./\./')
    tail -n +2 "$remaining_files" > $remaining_files.temp
    mv $remaining_files.temp $remaining_files
    # remove empty lines since they produce false positives
    sed '/^$/d' $fileA > tempA

    echo Comparing $fileA with other files...

    while read string; do
        fileB=$(echo $string | sed 's/.[^.]*\./\./')
        sed '/^$/d' $fileB > tempB
        A_len=$(cat tempA | wc -l)
        B_len=$(cat tempB | wc -l)

        differences=$(sdiff -B -s tempA tempB | wc -l)
        common=$(expr $A_len - $differences)

        percentage=$(echo "100 * $common / $B_len" | bc)
        if [[ $percentage -gt 15 ]]; then
            echo "  $percentage% duplication in" \
                 "$(echo $fileB | sed 's|\./||')"
        fi
    done < "$remaining_files"
    echo " "
done < "$all_files"

rm tempA
rm tempB
rm $all_files
rm $remaining_files
