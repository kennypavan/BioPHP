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
//prints 11
```

### Count Point Mutations Between Two Sequences
```php
$BioPHP = new BioPHP('CTGATGATGGGAGGAAATTTCA','CTGATGATGCGAGGGAATATCG');
echo $BioPHP->countPointMutations();
//prints 4
```
