# BioPHP - PHP Bioinformatics class
BioPHP is an open source project that implements some light tools for manipulating genomic data.

## Simple Usage:

```php
$BioPHP = new BioPHP('ATGAAa');
$BioPHP->reverseSequence();
$BioPHP->complementDnaSequence();
echo $BioPHP->sequenceA;
```
