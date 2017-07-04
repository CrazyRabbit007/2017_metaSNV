#! /bin/bash
############################################
#  metaSNP Step III:   `Post-Processing    #
############################################
#
# This code is part of the metagenomic SNP calling pipeline (metaSNP)
# Copyright (c) 2016 Robin Munch
# Licenced under the GNU General Public License (see LICENSE) 
#
# Helper script for metaSNP distance computation


display_usage() {
	echo >&2 ""
	echo >&2 "  Usage: $0 input_dir/ output/"
	echo >&2 ""
	echo >&2 "Parameters:"
	echo >&2 "	Required:"
	echo >&2 "	  input_dir/	DIR	directory with the filtered frequency tables (metaSNP_filtering.py)."
	echo >&2 "	  output_dir/	DIR	write output to DIR."	
	echo >&2 ""
	echo >&2 "Note: Expecting 'metaSNP_filtering.py' to be completed."
	echo >&2 ""
}

required_parameter() {
        echo >&2 ""
        echo >&2 "ERROR: '$1' is a required parameter"
        display_usage
        exit 1
}

dir_missing() {
        echo >&2 ""
        echo >&2 "ERROR: '$1' no such file or directory"
        display_usage
        exit 1
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
IN_DIR="$1"
OUT_DIR="$2"
R="$3"

## Test
[ -n "$IN_DIR" ] || required_parameter "input_dir/"
[ -n "$OUT_DIR" ] || required_parameter "output_dir/"

[ -d "$IN_DIR" ] || dir_missing "$IN_DIR"
[ -d "$OUT_DIR" ] || dir_missing "$OUT_DIR"

## Commandlines:
for i in $(ls $IN_DIR"/"*.freq); do echo "$R --no-save --file=$DIR/src/computeDistances.R --args $i $OUT_DIR"; done

exit
