#!/usr/bin/env bash

for file in $(find . -name 'CMakeLists.txt' && find cmake/ -name '*.cmake')
do
    cmake-format $file > .deleteme.txt
    diff -u --color=always $file .deleteme.txt
    diff -u $file .deleteme.txt >> .diff_deleteme.txt
done

if [ -s .deleteme.txt ]
then
    rm .deleteme.txt
fi

if [ -s .diff_deleteme.txt ]
then
    rm .diff_deleteme.txt
    exit 1
fi
