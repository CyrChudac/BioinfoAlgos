# BioinfoAlgos

This toolbox is capable to handle some fundamental bioinformatic tasks using the biopython library.

## How to run it?

Download the tool either as zip and extract it, or with the

      git init
      git pull https://github.com/CyrChudac/BioinfoAlgos
      
commands. Then continue with

      python -m venv yourEnvirName
     
here you can name your virtual environment as you want. Now you need to activate the virtual environment.
On Linux and MacOs simply run the *yourEnvirName/bin/activate* bash script,
on Windows run the *yourEnvirName\Scripts\activate.bat* file via command line.

Here you need to intall the required libraries using the

      pip install -r requirements.txt

command. Now simply run the python *main.py* script as you normally would:

      python main.py <and here the passed arguments>

## Structure

The only subdirectory in the repository is the directory [*testInputs*](https://github.com/CyrChudac/BioinfoAlgos/tree/master/testInputs)
where all of the testing files are stored,
there is also a file with arguments called [*testingInputs.txt*](https://github.com/CyrChudac/BioinfoAlgos/blob/master/testInputs/testingInputs.txt)
which is ready for Windows. If you wish to run it on Linux simply use the 

      cat ./testInputs/testingInputs.txt | tr "\\" "/" > ./testInputs/newName

command and then you have your file ready to test with the unix file separators.

## Architecture

Every functionality the software provides is encapsulated into a class of suitable name. 
All of those classes derive from a base class called Task in the [*Task.py*](https://github.com/CyrChudac/BioinfoAlgos/blob/master/Task.py) file.
They all return a Result class that is also defined in the *Task.py* file. Result remembers the value and also delegate to the function, that displays it.
That way I can easily call Tasks and just get their value without printing the result. Each task also has it's own help method, 
and so the help page is obtained by iterating through the whole dictionary and when a parameter is wrong deeper in the parameter structure, you
only see the help related to it and not the whole help page.


In the main.py the desired Task is found via dictionary with first arguments as keys. If the argument doesnt match any, the help task is assigned and help is printed.

## Functionalities

seen also with the -h command

### Fasta parsing

This functionality let's you parse the fasta files and see the sequences, their lengthses, their descriptions and subsequences.

### Hamming distance

Well... This functionality computes the hamming distance given two sequences (they have to have the same length).

### Edit distance

This functionality computes the edit distance of two sequences via dynamic programming. if specified, it can also display it's dynamic matrix and the final alignment. 
The final alignment is not computed (backtrack), unless -a or -ma is used.

### PDB parsing

This functionality parses the PDB files. It can show residues, atoms or both. It can count quantities of those. It can also count the width of structure = 
the disance between the far most atoms. It can find all atoms or residues around a ligand of given name. And you may also just display the 
list of ligands present in the structure.

### MSA parsing

This functionality lets you parse the MSA files created in the CLUSTAL format. It can simply display the alignment, show sequence given it's name or position,
show column given it's position or compute sum of pairs score for the whole alignement or given column. For the sum of pairs score a scoring matrix may be provided.
It's description is to be seen below.

### Position conservation

This functionality provides the n best scoring columns of an alignment provided a scoring matrix (it's description is below) and, obviously, the alignment itself. It can also determine whether are
residues around given position are conservated enough.

### MSA and structure combination

This functionality answers the question, if the residues around an active site of a ligand are conservated. 


## Scoring Matrix format
The scoring matrix is provided as a file in following format:

      default_missmatch|defaultMatch
      \tN\tA\tM\tE\tS\t
      N\t0\t1\t0.71\t0.12\t0.32\t
      A\t0.8\t-0.01\t0.3\t0.6\t0.198\t
          .
          .
          .

So the first line has default missmatch value, vertical bar and default match value.
the second line have names of residues separated by tab. This line starts and ends with a tab too.
All of the following lines start with a residue name followed by row of values, all separated by tabs.
