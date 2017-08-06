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
		return $this;

	}


	public function reverseSequence()
 	{

 		$this->sequenceA = strrev($this->sequenceA);
		return $this;

	}

	public function complementDnaSequence()
	{

		$this->sequenceA = str_replace("A", "t", $this->sequenceA);
		$this->sequenceA = str_replace("T", "a", $this->sequenceA);
		$this->sequenceA = str_replace("G", "c", $this->sequenceA);
		$this->sequenceA = str_replace("C", "g", $this->sequenceA);
		$this->normalizeSequence();
		return $this;

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
?>