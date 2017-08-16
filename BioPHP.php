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

	public $codonToAminos = ["ATT" => "I", "ATC" => "I", "ATA" =>"I", "CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L", "TTA" => "L", "TTG" => "L", "GTT" => "V", "GTC" => "V", "GTA" =>"V", "GTG" =>"V","TTT" => "F", "TTC" => "F","ATG" => "M","TGT" => "C", "TGC" => "C", "GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A", "GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G", "CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P", "ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T", "TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S", "AGT" => "S", "AGC" => "S", "TAT" => "Y", "TAC" => "Y", "TGG" => "W", "CAA" => "Q", "CAG" => "Q", "AAT" => "N", "AAC" => "N", "CAT" => "H", "CAC" => "H", "GAA" => "E", "GAG" => "E", "GAT" => "D", "GAC" => "D", "AAA" => "K", "AAG" => "K", "CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R", "AGA" => "R", "AGG" => "R", "TAA" => "*", "TAG" => "*", "TGA" => "*" ]; // * equals stop codon

	public $monisotopicAminoMass = ["A" => 71.03711, "C" => 103.00919, "D" => 115.02694, "E" => 129.04259, "F" => 147.06841, "G" => 57.02146, "H" => 137.05891, "I" => 113.08406, "K" => 128.09496, "L" => 113.08406, "M" => 131.04049, "N" => 114.04293, "P" => 97.05276, "Q" => 128.05858, "R" => 156.10111, "S" => 87.03203, "T" => 101.04768, "V" => 99.06841, "W" => 186.07931, "Y" => 163.06333];


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
		$results = [];

		for($i=0; $i<=$sLen; $i++)
		{

			if(substr($this->sequenceB, $i, $tLen) == $this->sequenceA){
				$results[] = $i+1;
			}

		}

		return implode(" ",$results);

	}


	public function getReadingFrames()
	{
		
		
		$frameOne = $this->sequenceA;
		$frameTwo = substr($this->sequenceA, 1);
		$frameThree = substr($this->sequenceA, 2);
		$readingFrames = [$frameOne,$frameTwo,$frameThree];

		return $readingFrames;

	}


	public function calcMonoIsotopicMass($proteinSequence)
	{
		
		$proteinLen = strlen($proteinSequence);
		$mass = 0;

		for($i=0; $i<=$proteinLen; $i++)
		{
			
			if( isset( $this->monisotopicAminoMass[substr($proteinSequence, $i, 1)] ) ){
				$mass = $mass + $this->monisotopicAminoMass[substr($proteinSequence, $i, 1)];
			}

		}		


		return $mass;

	}


	//read fasta as a string and return an array.
	public function readFasta($fastaStr)
	{
		
		$fastaLines = explode('>', $fastaStr);
		$fastaArray = [];

		foreach ($fastaLines as $fastaLine) 
		{

			$singleLines = preg_split('/$\R?^/m', $fastaLine);

			$sequence = "";

			for ($i=1;$i<count($singleLines);$i++) 
			{
			
				$sequence .= str_replace(array("\r", "\n"), '',$singleLines[$i]);
			
			}
			
			$fastaArray[] = ["name" => $singleLines[0], "sequence" => $sequence];
			
		}

		return $fastaArray;

	}



	public function mostLikelyCommonAncestor($sequencesArray)
	{
		
		$countNucleotides = [];

		foreach ($sequencesArray as $key => $value) {
			
			$sequence = str_split($value['sequence']);

			for ($i=0; $i<count($sequence);$i++){

				if(!array_key_exists($i, $countNucleotides)){
					$countNucleotides[$i]["A"] = 0;
					$countNucleotides[$i]["T"] = 0;
					$countNucleotides[$i]["G"] = 0;
					$countNucleotides[$i]["C"] = 0;
				}

				if($sequence[$i] == "A"){
					$countNucleotides[$i]["A"]++;
				} elseif($sequence[$i] == "T"){
					$countNucleotides[$i]["T"]++;
				} elseif($sequence[$i] == "G"){
					$countNucleotides[$i]["G"]++;
				} elseif($sequence[$i] == "C"){
					$countNucleotides[$i]["C"]++;
				}				

			}

		}

		$mostLikelyCommonAncestorSequence = "";
		
		for($i=0;$i<count($countNucleotides);$i++)
		{

			$mostLikelyCommonAncestorSequence .= array_search(max($countNucleotides[$i]), $countNucleotides[$i]);

		}

		return $mostLikelyCommonAncestorSequence;

	}


	public function getUniprotFastaByID($uniprotID)
	{
		
		
		$ch = curl_init(); 
		curl_setopt($ch, CURLOPT_URL, "http://www.uniprot.org/uniprot/".$uniprotID.".fasta"); 
		curl_setopt($ch, CURLOPT_FOLLOWLOCATION, true);
		curl_setopt($ch, CURLOPT_RETURNTRANSFER, 1); 
		$output = curl_exec($ch); 
		curl_close($ch);  		

		return $output;

	}

	//create and array of all matchable amino acids at each position.
	public function varyingFormsGeneration($varyingSubSequence)
	{
		
		$varyingSubSequence = str_split($varyingSubSequence);
		$squareBrace = false;
		$curlyBrace = false;
		$returnedVaryingSubSequence = [];

		$inc=0;

		for ($i=0; $i < count($varyingSubSequence); $i++) { 
			
			if($varyingSubSequence[$i] == "]"){

				$squareBrace = false;

			} elseif($varyingSubSequence[$i] == "[" || $squareBrace == true){

				if($squareBrace == true){

					$returnedVaryingSubSequence[$inc][] = $varyingSubSequence[$i];

				}

				$squareBrace = true;

				continue;

			} elseif($varyingSubSequence[$i] == "}"){
				
				$curlyBrace = false;

			} elseif($varyingSubSequence[$i] == "{"  || $curlyBrace == true){

				if($curlyBrace == true){

					$returnedVaryingSubSequence[$inc][] = "!".$varyingSubSequence[$i];

				}

				$curlyBrace = true;

				continue;

			} else {

				$returnedVaryingSubSequence[$inc] = $varyingSubSequence[$i];
				
			}

			$inc++;
		}

		return $returnedVaryingSubSequence;

	}


	public function findMotifProtein($varyingSubSequence,$proteinID)
	{

		// find the variations in subsequence
		$varyingSubSequences = $this->varyingFormsGeneration($varyingSubSequence);

		// get sequence from uniprot
		$uniprotFasta =  $this->getUniprotFastaByID($proteinID);
		$fastaArray = $this->readFasta($uniprotFasta);
		$sequence = $fastaArray[1]['sequence'];

		// get sequence lengths and declare results
		$tLen = count($varyingSubSequences);
		$sLen = strlen($sequence);
		$results = [];

		// search for motifs in the sequence
		for($i=0; $i<=($sLen-$tLen); $i++)
		{

			$matches = 0;

			for($b=0; $b<$tLen; $b++)
			{

				if(is_array($varyingSubSequences[$b])){

					foreach ($varyingSubSequences[$b] as $singleValue) 					
					{

						if($singleValue[0] == "!"){

							if($singleValue[1] != $sequence[$i+$b]){

								$matches++;

							}

						} else {

							if($singleValue == $sequence[$i+$b]){

								$matches++;

							} 

						}

					}

				} else {

					if($varyingSubSequences[$b] == $sequence[$i+$b]){

						$matches++;

					} 

				}

				if($matches == $tLen){

					$results[] = ($i+1);
				
				}

			}

		}

		return $results;

	}

}


//--------- Sample Usage Below ------------//

/*
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

//Sample Usage - Get all reading frames (in one direction)
$BioPHP = new BioPHP('GTATATCTATATGGCCATAT');
print_r( $BioPHP->getReadingFrames() );
echo "\n";

//Sample Usage - Calculate monoisotopic mass
$BioPHP = new BioPHP('CTGATGATGGGAGGAAATTTCAGA');
$proteinSequence = $BioPHP->translateDna()."\n";
echo $BioPHP->calcMonoIsotopicMass($proteinSequence)."\n\n";


//Sample Usage - Find most common likely ancestor
$fasta = "
>Sequence 1
ATCCAGCT
>Sequence 2
GGGCAACT
>Sequence 3
ATGGATCT
";
$BioPHP = new BioPHP();
$fastaArray = $BioPHP->readFasta($fasta);
echo $BioPHP->mostLikelyCommonAncestor($fastaArray)."\n";


//Sample Usage - Get a fasta result from Uniprot and calculate its Isotpoic Mass
$BioPHP = new BioPHP();
$uniprotFasta =  $BioPHP->getUniprotFastaByID("B5ZC00");
$fastaArray = $BioPHP->readFasta($uniprotFasta);
echo $BioPHP->calcMonoIsotopicMass($fastaArray[1]['sequence'])."\n\n";


//Sample Usage - Get a protein fasta result from Uniprot and find protein motif with varying sequence search.
// Varying sequence - [XY] means "either X or Y" and {X} means "any amino acid except X."  N-glycosylation motif is written as N{P}[ST]{P}.		
$BioPHP = new BioPHP();
$results = $BioPHP->findMotifProtein("N{P}[ST]{P}","B5ZC00");
print_r($results);
*/
?>