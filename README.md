# BioPHP - PHP Bioinformatics class
BioPHP is an easy to use, open source project. BioPHP implements a selection of simple tools for manipulating genomic data. We're not doing any heavy lifting here (this is PHP). This class is built for basic RNA, DNA, and Protein manipulation. 

## Simple Usage:

### Find Reverse Complement
```php
$BioPHP = new BioPHP('ATGAAA');
$BioPHP->reverseSequence();
$BioPHP->complementDnaSequence();
echo $BioPHP->sequenceA;
//prints TTTCAT
```

### Calculate GC Content
```php
$BioPHP = new BioPHP('ATGAAAGCATC');
echo $BioPHP->gcContent();
//prints 36.3636
```

### Count Point Mutations Between Two Sequences
```php
$BioPHP = new BioPHP('CTGATGATGGGAGGAAATTTCA','CTGATGATGCGAGGGAATATCG');
echo $BioPHP->countPointMutations();
//prints 4
```

### Translate DNA Sequence to Amino Acid Sequence
```php
$BioPHP = new BioPHP('CTGATGATGGGAGGAAATTTCAGA');
echo $BioPHP->translateDna();
//prints LMMGGNFR
```

### Finding a Motif in DNA
```php
$BioPHP = new BioPHP('ATAT', 'GTATATCTATATGGCCATAT');
echo $BioPHP->findMotifDNA();
//prints 3 9 17
```

### Get all three reading frames
```php
$BioPHP = new BioPHP('GTATATCTATATGGCCATAT');
print_r( $BioPHP->getReadingFrames() );

/*
* returns array containing...
Array
(
    [0] => GTATATCTATATGGCCATAT
    [1] => TATATCTATATGGCCATAT
    [2] => ATATCTATATGGCCATAT
)
*/
```
