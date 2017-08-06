<?php

/**
 * BioPHP
 *
 * @category   BioInformatics
 * @package    BioPHP
 * @author     Kenny Pavan  <kenpavan@gmail.com>
 * @license    MIT
 */

class BioPHP {


	public $sequenceA;

	public $sequenceB;
	
	public $codonToAminos = array("ATT" => "I", "ATC" => "I", "ATA" =>"I","CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L", "TTA" => "L", "TTG" => "L", "GTT" => "V", "GTC" => "V", "GTA" =>"V", "GTG" =>"V","TTT" => "F", "TTC" => "F","ATG" => "M","TGT" => "C", "TGC" => "C", "GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A", "GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G", "CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P", "ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T", "TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S", "AGT" => "S", "AGC" => "S", "TAT" => "Y", "TAC" => "Y", "TGG" => "W", "CAA" => "Q", "CAG" => "Q", "AAT" => "N", "AAC" => "N", "CAT" => "H", "CAC" => "H", "GAA" => "E", "GAG" => "E", "GAT" => "D", "GAC" => "D", "AAA" => "K", "AAG" => "K", "CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R", "AGA" => "R", "AGG" => "R", "TAA" => "*", "TAG" => "*", "TGA" => "*" ); // * equals stop codon


	public function __construct($sequenceA=false, $sequenceB=false)
	{
		$this->sequenceA = $sequenceA;
		$this->sequenceB = $sequenceB;
		$this->normalizeSequence();
	}


	public function normalizeSequence()
	{

 		$this->sequenceA = strtoupper($this->sequenceA);
 		$this->sequenceB = strtoupper($this->sequenceB);
		return $this->sequenceA;

	}


	public function reverseSequence()
 	{

 		$this->sequenceA = strrev($this->sequenceA);
		return $this->sequenceA;

	}


	public function complementDnaSequence()
	{

		$this->sequenceA = str_replace("A", "t", $this->sequenceA);
		$this->sequenceA = str_replace("T", "a", $this->sequenceA);
		$this->sequenceA = str_replace("G", "c", $this->sequenceA);
		$this->sequenceA = str_replace("C", "g", $this->sequenceA);
		$this->normalizeSequence();
		return $this->sequenceA;

	}


	public function countNucleotides()
	{

		return strlen($this->sequenceA);

	}


	public function gcContent($precision = 2)
	{


		$g = substr_count($this->sequenceA, 'G');
		$c = substr_count($this->sequenceA, 'C');

		return number_format( ( ($g+$c) / strlen($this->sequenceA) ) * 100 , $precision);

	}


	public function convertRnaToDna()
	{

		$this->sequenceA = str_replace("U","T",$this->sequenceA);
		return $this->sequenceA;

	}


	public function countPointMutations()
	{

		$totalMutations = 0;

		if(strlen($this->sequenceA) >= strlen($this->sequenceB)){

			$longestSequenceLength = strlen($this->sequenceA);

		} else {
			
			$longestSequenceLength = strlen($this->sequenceB);

		}
		

		for($i=0; $i < $longestSequenceLength; $i++)
		{

			if (isset($this->sequenceA[$i]) && isset($this->sequenceB[$i]) ) {

				if($this->sequenceA[$i] != $this->sequenceB[$i]) {
					$totalMutations++;
				}				

			} else {

				$totalMutations++;

			}

		
		}

		return $totalMutations;

	}


	public function translateDna($offset = 0){

		$sequence = substr($this->sequenceA, $offset);  //offset reading frame for future use.
		$sequenceCodons = str_split($sequence, 3);

		$proteinSequence = "";

		foreach ($sequenceCodons as $sequenceCodon) 
		{
		  
			if(isset($this->codonToAminos[$sequenceCodon])){

				$proteinSequence .= $this->codonToAminos[$sequenceCodon]; 

			} else {

				$proteinSequence .= "-"; //unknown

			}

		}

		return $proteinSequence;
	}

	//sequence A is a substring of sequence B
	public function findMotifDNA()
	{

		$tLen = strlen($this->sequenceA);
		$sLen = strlen($this->sequenceB);

		for($i=0; $i<=$sLen; $i++)
		{
			
			if(substr($this->sequenceB, $i, $tLen) == $this->sequenceA){
				$results[] = $i+1;
			}

		}

		return implode(" ",$results);

	}

}


//Sample Usage - Find Reverse Complement
$BioPHP = new BioPHP('ATGAAAGCATC');
$BioPHP->reverseSequence();
$BioPHP->complementDnaSequence();
echo $BioPHP->sequenceA."\n";


//Sample Usage - Get Nucleotide Count
$BioPHP = new BioPHP('ATGAAAGCATC');
echo $BioPHP->countNucleotides()."\n";


//Sample Usage - Get GC Content Percentage
$BioPHP = new BioPHP('ATGAAAGCATC');
echo $BioPHP->gcContent(4)."\n";


//Sample Usage - Count point Mutations between two sequences
$BioPHP = new BioPHP('CTGATGATGGGAGGAAATTTCA','CTGATGATGCGAGGGAATATCG');
echo $BioPHP->countPointMutations()."\n";


//Sample Usage - Convert RNA to DNA
$BioPHP = new BioPHP('ACGCGAUUGCGAUCGAUGCACGCU');
echo $BioPHP->convertRnaToDna()."\n";


//Sample Usage - Translate sequence to amino acid
$BioPHP = new BioPHP('CTGATGATGGGAGGAAATTTCAGA');
echo $BioPHP->translateDna()."\n";


//Sample Usage - Finding a Motif in DNA
$BioPHP = new BioPHP('ATAT', 'GTATATCTATATGGCCATAT');
echo $BioPHP->findMotifDNA()."\n";


?>