# BioPHP - PHP Bioinformatics class
BioPHP is an easy to use, open source project. BioPHP implements a selection of simple tools for manipulating genomic data. We're not doing any heavy lifting here (this is PHP). This class is built for basic RNA, DNA, and Protein manipulation. 

## Simple Usage:

```php
$BioPHP = new BioPHP('ATGAAa');
$BioPHP->reverseSequence();
$BioPHP->complementDnaSequence();
echo $BioPHP->sequenceA;
// prints TTTCAT
```
