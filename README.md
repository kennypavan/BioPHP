# BioPHP - PHP Bioinformatics class
BioPHP is an easy to use, open source project, which implements light tools for manipulating genomic data. We're not doing any heavy lifting here. This PHP class is built for basic RNA, DNA, and Protein manipulation. 

## Simple Usage:

```php
$BioPHP = new BioPHP('ATGAAa');
$BioPHP->reverseSequence();
$BioPHP->complementDnaSequence();
echo $BioPHP->sequenceA;
// prints TTTCAT
```
