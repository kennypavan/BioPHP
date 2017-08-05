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
	public $sequenceC;


	public function __construct($sequenceA=false, $sequenceB=false, $sequenceC=false)
	{
		$this->sequenceA = $sequenceA;
		$this->normalizeSequence();
	}


	public function normalizeSequence()
	{

 		$this->sequenceA = strtoupper($this->sequenceA);
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

}

$BioPHP = new BioPHP('ATGAAa');
$BioPHP->reverseSequence();
$BioPHP->complementDnaSequence();
echo $BioPHP->sequenceA;
?>