#!/bin/bash
FASTA="$1"
FILENAME=$(basename "$FASTA")

echo -e "Filename\tA\tT\tG\tC\tN\tOther"

if [[ "$FASTA" == *.gz ]]; then
  zcat "$FASTA" \
  | awk -v file="$FILENAME" 'BEGIN{A=T=G=C=N=other=0}
    /^[>]/ {next}
    {
      line=toupper($0)
      A+=gsub(/A/,"",line)
      T+=gsub(/T/,"",line)
      G+=gsub(/G/,"",line)
      C+=gsub(/C/,"",line)
      N+=gsub(/N/,"",line)
      other += length(line) - (A+T+G+C+N - other) # careful accumulation
    }
    END{printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\n",file,A,T,G,C,N,other}'
else
  awk -v file="$FILENAME" 'BEGIN{A=T=G=C=N=other=0}
    /^[>]/ {next}
    {
      line=toupper($0)
      A+=gsub(/A/,"",line)
      T+=gsub(/T/,"",line)
      G+=gsub(/G/,"",line)
      C+=gsub(/C/,"",line)
      N+=gsub(/N/,"",line)
      other += length(line) - ( (gsub(/A/,"",line)) + (gsub(/T/,"",line)) + (gsub(/G/,"",line)) + (gsub(/C/,"",line)) + (gsub(/N/,"",line)) )
    }
    END{printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\n",file,A,T,G,C,N,other}' "$FASTA"
fi
