#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd "$SCRIPT_DIR" || exit
rm -f main.pdf 2> /dev/null
rm -f ./src/chapters/*.aux
mv ./src/main.pdf . 2> /dev/null
cd ./src/ || exit
GLOBIGNORE=*.tex:*.bib
rm -v * 2> /dev/null