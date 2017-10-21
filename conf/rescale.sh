#!/bin/sh

if [ $# -ne 2 ]
then
    echo "Usage: [file] [geometry]\nwhere file is a pgm file and geometry is like 200x40" >&2
    exit
fi

convert $1 -filter point -resize $2 -compress none ${1%.pgm}_$2.pgm
