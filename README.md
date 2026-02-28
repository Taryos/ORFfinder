ORF Finder

Description :

Script Python permettant d'identifier les cadres de lecture ouverts (ORF) dans des séquences d'ADN au format FASTA.

Utilisation :

Lancez le script via le terminal avec la commande suivante :

```Bash 
python3 ORFfinder.py -i input.fa -o output.fa [options]
```

Options principales :

    -m : Longueur minimale de l'ORF (défaut : 100).

    -t : Active la traduction en acides aminés.

    -r : Analyse également le brin reverse-complément.

    -s : Utilise uniquement ATG comme codon Start (Start strict).

   

