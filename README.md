# BioPHP - PHP Bioinformatics class
BioPHP is an easy to use, open source project. BioPHP implements a selection of simple tools for manipulating genomic data. We're not doing any heavy lifting here (this is PHP). This class is built for basic RNA, DNA, and Protein manipulation. 

## Simple Usage:

### Find Reverse Complement
```php
$BioPHP = new BioPHP('ATGAAA');
$BioPHP->reverseSequence(); //reverse sequence
$BioPHP->complementDnaSequence(); //get the reversed complement
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

### Calculate Monoisotopic Mass
```php
$BioPHP = new BioPHP('CTGATGATGGGAGGAAATTTCAGA');
$proteinSequence = $BioPHP->translateDna()."\n"; //translate sequence
echo $BioPHP->calcMonoIsotopicMass($proteinSequence)."\n"; //calculate mass
//prints 906.42041
```

### Finding a Motif in DNA
```php
$BioPHP = new BioPHP('ATAT', 'GTATATCTATATGGCCATAT');
echo $BioPHP->findMotifDNA();
//prints 3 9 17
```

### Get Reading Frames
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

//Protip: To get all 6 reading frames. Use the reverse and complement methods, then pass the result to getReadingFrames()
```


### Find most common likely ancestor
```php
$fastaSequence = "
>Sequence 1
ATCCAGCT
>Sequence 2
GGGCAACT
>Sequence 3
ATGGATCT
";

$BioPHP = new BioPHP();
$fastaArray = $BioPHP->readFasta($fastaSequence); //read and parse the sequences
echo $BioPHP->mostLikelyCommonAncestor($fastaArray)."\n";

//prints ATGCAACT
```
