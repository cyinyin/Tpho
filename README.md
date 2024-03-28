# Tpho :Predicting optimal catalytic temperature of phosphatases by removing conserved amino acids in phosphatases using machine-learning

Using machine learning to reveal that non-conserved amino acids affect the optimal catalytic temperature of enzymes to a greater extent than conserved amino acids using phosphatase as an example. 

Phosphatase is a hydrolase that belongs to the halogen dehalogenase superfamily, he is an in vitro polymerase that catalyzes the last step, irreversible, thermodynamically prone reaction that can pull the whole system towards the production of a product. In this work, machine learning (ML) is used to associate the amino acid sequence of the phosphatase with the optimal catalytic temperature by identifying key sequence features associated with the optimal catalytic temperature of the phosphatase utilized by the ML algorithm.The strategy used in this work may be applicable to reveal sequence-function relationships in other protein families.   
Tpho is a Pyhton toolkit for predicting the optimal catalytic temperature of phosphatases.Tpho uses a regressor to predict the optimal catalytic temperature of phosphatases.
## Requirements

- Python(>=3)
    
Python modules（version used in this work）    

- pandas(1.3.2)
- numpy(1.21.5)
- scipy(1.7.3)
- biopython(1.79)
- scikit-learn(0.24.1)
- matplotlib(3.4.2)
- pydpi(1.0)
- keras(2.3.1)
- requests(2.26.0)
 
## Project structure
```bash
root
|-- data        // Storing data
|-- get_dataset       // Access to data information
|   |-- Spider_Phosphatase_EC_Numbre.py     // Climbe phosphatase EC number
|   |-- Spider_Phosphatase_Information_Html.py     // Crawle phosphatase information page
|   |-- Parsing_Phosphatase_Information_Html.py  // Parse phosphatase information page 
|-- predicate_phosphatase_temperature // Predicte phosphatase temperature
|   |-- AAComposition.py     // Calculate amino acid composition
|   |-- AAIndex.py     // obtain the properties of amino acids or their pairs from the aaindex database
|   |-- Autocorrelation.py  //  Calculate the Autocorrelation descriptors based different properties of AADs 
|   |-- ConjointTriad.py     //  Calculate the conjoint triad features
|   |-- CTD.py     // Calculate composition, transition, distribution
|   |-- GetSubSeq.py  // Split the total protein into a set of segments around specific amino 
|   |-- PseudoAAC.py  // Calculate pseudo amino acid composition 
|   |-- PyPro.py     // A class used for computing different types of protein descriptors!
|   |-- QuasiSequenceOrder.py     // Calculate quasi-sequence order
|   |-- resreg.py  // Regression resampling 
|   |-- Tpho.py  // Predicte the optimal catalytic temperature of phosphatase 
|-- processing_data    // Data processing
|   |-- clustal-omega-1.2.2-win64     // Multiple Sequence Comparison Tool
|   |-- CompleteInformationSequence.py     // Obtain the original sequence information corresponding to the de-conserved amino acid sequence
|   |-- MultipleSequenceAlignment.py     // Multiple Sequence Comparison 
|   |-- ProcessingComparedSequence.py  // De-conserved amino acids in the original sequence 
|-- result_validation // Further validation of prediction performance
|   |-- collationsequence.py     // Organize validation data information
|   |-- resultvalidation.py     // Predicting results using application model
|-- test // Test model
|   |-- test_data     // Test data information
|   |-- test_model.py     // Test  application model
|-- unit         // Functional unit
|   |-- path.py  // Path information
```
## Python Scripts
Main scripts (in chronological order as used in the study)  
- MultipleSequenceAlignment.py: Multiple Sequence Comparison with clustalo.exe.
- ProcessingComparedSequence.py: Removal of conserved amino acids from the original sequence.
- Tpho.py:Predicting the optimal catalytic temperature of phosphatases using machine learning.
## Usage    
You need to first perform multiple sequence comparisons, create a new FASTA folder, place the predicted phosphatase and training set data in this folder, and use the online tool clustalw（ https://www.genome.jp/tools-bin/clustalw ）Perform sequence similarity comparison and download the DND file. Subsequently, the result DND file was visualized using MEGA11 software to construct an evolutionary tree. Determine the position of the phosphatase Uniprot id in the evolutionary tree and determine the phosphatase sequence EC number. Then, put the comparison phosphatase information into the original_ Complete_ Information_ In sequence. tsv, in addition, place this sequence and phosphatases under the same EC number in fasta format in the ComparePathatease file, run the MultipleSequenceAlignment. py file to obtain the sequence_ Alignment_ Conservative amino acid index information in the results. tsv file, then place the information in the sequence_ Alignment_ Results_ Intermediate_ File.tsv, run the ProcessingComparedSequence.py file, remove the conserved amino acids from the sequence, and place the results in test_data.tsv file. Put the phosphatase information with conserved amino acids removed into the test_data.fasta file. run the test_model.py file.
## Examples:
- ```text
  
  A0A2W4KKP5		3.1.3.102
  B9L0J9     3.1.3.64
    
  sequence_alignment_results.tsv
  A0A2W4KKP5		[9, 11, 15, 50, 87, 117, 129, 134, 137, 164, 165, 172, 181, 184, 201]		['L', 'D', 'T', 'W', 'D', 'L', 'N', 'R', 'I', 'K', 'P', 'E', 'P', 'C', 'G']
  B9L0J9     [17, 36, 126, 169, 219, 292, 293]     ['E', 'G', 'P', 'R', 'L', 'R', 'W']
        
  from predicate_phosphatase_temperature import Tpho
  fasta_file = 'test/test_data/experimental_data.fasta'
  pred_data = Tpho.pred_phosphatase_topt(fasta_file)
  print（pred_data）   
        uniprot id topt_pred
    0   A0A2W4KKP5      70.0
    1   B9L0J9      80.0
  ```
