#!/bin/bash
id=$$
if [ ! -z "$1" ] ; then
  id=$1
fi
id="$id"`date +%s`
aaseq=`cat -`
echo $aaseq > webserv$id.fasta
./asaquick webserv$id.fasta
echo "ASA and its ERR are in Angstrom squared"
echo 1 | awk '{printf("%-6s%2s %6s %6s %8s %8s\n", "ID", "AA", "ASA", "ERR", "RASA", "ERR")}'
paste asaq.webserv$id.fasta/asaq.pred asaq.webserv$id.fasta/rasaq.pred | awk '{printf("%-6s %1s %6.0f %6.0f %8.3f %8.3f\n",$1,$2,$3,$4,$7,$8)}'

exit


