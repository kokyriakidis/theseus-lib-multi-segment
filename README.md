# Theseus-lib

**TABLE OF CONTENTS:**

* [Introduction](#introduction)
    * [What is Pericles?](#what_is_pericles)
    * [Getting started](#getting_started)
* [Using Theseus](#using_theseus)
    * [Consensus Generation](#consensus_usage)
* [Tools](#theseus_tools)
* [Heuristics](#theseus_heuristics)
* [Reporting Bugs and Feature Request](#theseus_bugs)
* [License](#theseus_licence)
* [Authors](#theseus_authors)

## <a name="introduction"></a> 1. Introduction

### <a name="what_is_pericles"></a> 1.1. What is Pericles?
Pericles isa specialized version of Theseus tailored for the use case of consensus generation via Partial Order Alignment (POA). Theseus is a fast, optimal and affine-gap Sequence-to-Graph aligner. It leverages the expected high similarity in the aligned data to accelerate computation and reduce the search space compared to other alternatives.

<p align = "center">
<img src = "img/Theseus_green.png" width="300px">
</p>

Theseus extends the proposal from the Wavefront Alignment algorithm (WFA), originally devised for pairwise sequence alignment, to the context of sequence-to-graph alignment.

### <a name="getting_started"></a> 1.2. Getting started
Git clone and compile the library, tools, and examples (by default, use `cmake` for the library and benchmark build).

```
git clone https://github.com/albertjimenezbl/theseus-lib
cd theseus-lib
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```


## <a name="using_theseus"></a> 2. Using Theseus

### <a name="consensus_usage"></a> 2.1. Consensus Generation

This example illustrates how to use Theseus as a Consensus generator. First, include the msa and general headers:
```
#include "theseus/alignment.h"
#include "theseus/heuristics.h"
#include "theseus/penalties.h"
#include "theseus/theseus_msa_aligner.h"
```

Then, create and configure a MSA aligner object. This object is defined by three parameters: a set of penalties, a heuristics object, and an initial sequence behaving as the starting point of the alignment. An example on how to set such parameters and create an MSA object is found in the next code snippet:
```
theseus::Penalties penalties(match, mismatch, gap_open, gap_extend);
theseus::Heuristics heuristics(use_lag_pruning, use_density_drop);
theseus::TheseusMSA aligner(penalties, heuristics, initial_sequence, is_ends_free);
```

Once this is done, we can start adding sequences to our POA graph using the align functionality, that returns an Alignment object with CIGAR, path, and score information. Each time that you call the aligner, you have to provide three parameters: 1) The sequence to be aligned to the poa graph, 2) A boolean variable indicating whether you want the alignment to go forward (source to sink) or reverse (sink to source), and 3) A boolean variable indicating whether you want the alignment to be end to end or having a free end at the "end" of the graph:
```
theseus::Alignment alignment_object = aligner.align(sequence, align_reverse, is_ends_free);
```

Each time a new sequence is added to the POA graph, the graph is updated with the newly found variation (all the insertions, deletions and mismatches of the resulting alignment object).

Finally, we can output the result in four different formats: a graph in .gfa format, a multiple sequence alignment, a consensus sequence, and a graph in .dot format:
```
aligner.print_as_gfa(output_file);             // Output the compacted POA graph in .gfa format
aligner.print_as_msa(output_file);             // Output as a Multiple Sequence Alignment
sequence = aligner.get_consensus_sequence();   // Find the consensus sequence of the alignment
aligner.print_as_dot(output_file);             // Output the compacted POA graph in .dot format
```


## <a name="theseus_tools"></a> 3.Tools

The Theseus library implements a minimal tools to use the Theseus algorithm on the consensus generation context. It is important to note that this tool is not production ready.

### <a name="consensus_tool"></a> 3.1. Consensus tool: pericles

This example illustrates how to use the **pericles** tool. This tool computes the MSA of the set of sequences in an given input *.fasta* file, allowing to add partial sequences, as long as they start on either end of a backbone sequence. The executable is located in the path */build/tools/pericles*:
```
cd build/tools/
```
**[IMPORTANT]**
The *.fasta* file containing sequences has a special structure. As all .fasta files, the data associated to each sequence has two parts: 1) A line starting with ">" containing metadata, and 2) the sequence itself, that appears on the next lines.

1) Constists of two elements: ">" and two boolean values (0 or 1) indicating 1) whether the associated read should be mapped in the canonical direction (from left to right and from source to sink) or in the reverse direction (from right to left and from sink to source), and 2) whether the added sequence has one free end on the graph or not.
2) Contains the sequence itself.

Examples aligning in the canonical direction using the end-to-end alignment and in the reversed direction using the ends free mode:
```
> 0 0
ACCT
> 1 1
TCCAT
```

Select the scoring scheme, set the input and output files and execute the tool. Each execution of *pericles* lets you select the following parameters:
```
Usage: pericles [OPTIONS]
                 Options:\n"
                   -m, --match <int>           The match penalty                                       [default=0]
                   -x, --mismatch <int>        The mismatch penalty                                    [default=2]
                   -o, --gapo <int>            The gap open penalty                                    [default=3]
                   -e, --gape <int>            The gap extension penalty                               [default=1]
                   -t, --output_type <int>     The output format of the multiple alignment             [default=0=MSA]
                                                0: MSA: Standard Multiple Sequence Alignment format,
                                                1: GFA: Output the resulting POA graph in GFA format,
                                                2: Consensus: Output the consensus sequence,
                                                3: Dot: Output in .dot format for visualization purposes.
                                                        Only tractable for small graphs
                   -f, --output <file>         Output file                                             [Required]
                   -s, --sequences <file>      Dataset file                                            [Required]"

                  Heuristics:\n"
                   -d  --density_heuristic     Activate the drop heuristic based on advancement density.
                   -l  --lag_pruning           Activate the pruning of diagonals lagging behind in the alignment.
```

An example of the execution of *pericles* is shown in the following piece of code
```
./pericles -m 0 -x 2 -o 3 -e 1 -t 0 -f output_file.out -s sequences.fasta
```


## <a name="theseus_heuristics"></a> 4. HEURISTICS
Theseus library implements some heuristic approaches that accelerate alignment at the expense of a limited loss in accuracy. In particular, Theseus implements 1) a **pruning heuristic** that discards diagonals that have fallen behind in the alignment, as long as the alignment has shown a significant advancement in the last scores, and 2) a **drop heuristic** that drops alignment when the advancement density (number of offsets advanced in the last scores) is very low.

You can activate these heuristics by initializing the Heuristics object with two separate boolean flags:
```
theseus::Heuristics heuristics(use_lag_pruning, use_density_drop);
```

Moreover, the three minimal tools provided in this alignment library allow you to activate these two heuristics from the command line, by adding the **-d** and **-l** flags:
```
./pericles -m 0 -x 2 -o 3 -e 1 -t 0 -f output_file.out -s sequences.fasta -d -l
```


## <a name="theseus_bugs"></a> 6. REPORTING BUGS AND FEATURE REQUEST

Feedback and bug reporting is highly appreciated. Please report any issue or suggestion on github or email to the main developer (albert.jimenez1@bsc.es). Don't hesitate to contact us if:
  - You experience any bug or crash.
  - You want to request a feature or have any suggestion.
  - Your application using the library is running slower than it should or you expected.
  - Need help integrating the library into your tool.


## <a name="theseus_licence"></a> 7. LICENSE

Theseus-lib is distributed under MIT licence.


## <a name="theseus_authors"></a> 8. AUTHORS

Albert Jimenez Blanco (albert.jimenez1@bsc.es) is the main developer and the person you should address your complaints.

Lorién López-Villellas has had major contributions both in the technical implementation of Theseus and the final structure of the library.

<!-- ## <a name="theseus_cite"></a> 7. CITATION

**Albert Jimenez-Blanco, Lorien Lopez-Villellas, Juan Carlos Moure, Miquel Moreto, Santiago Marco-Sola**. ["Theseus: Fast and Optimal Affine-Gap Sequence-to-Graph Alignment"](). Bioinformatics, 2026. -->

