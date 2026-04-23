# Theseus-lib

**TABLE OF CONTENTS:**

* [Introduction](#introduction)
    * [What is Theseus?](#what_is_theseus)
    * [Getting started](#getting_started)
* [Using Theseus](#using_theseus)
    * [Multiple Sequence Alignment (MSA)](#msa_usage)
    * [Aligning a sequence to a graph](#alignment_usage)
* [Tools](#theseus_tools)
* [Reporting Bugs and Feature Request](#theseus_bugs)
* [License](#theseus_licence)
* [Authors](#theseus_authors)

## <a name="introduction"></a> 1. Introduction

### <a name="what_is_theseus"></a> 1.1. What is Theseus?
Theseus is a fast, optimal and affine-gap Sequence-to-Graph aligner. It leverages the expected high similarity in the aligned data to accelerate computation and reduce the search space compared to other alternatives. Theseus is a general purpose aligner, providing two main functionalities:
1. **Multiple Sequence Alignment (MSA):** Theseus performs MSA of a set of N sequences using the Partial Order Alignment (POA) approach. That is, it progressively builds a partial order graph representing an MSA of a set of given sequences, adding one more sequence to the graph at each iteration.
2. **Alignment to a reference graph:** Alternatively, Theseus can align one or several sequences to a reference graph. The user has to provide an initial position for the alignment to start.

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

### <a name="msa_usage"></a> 2.1. Multiple Sequence Alignment (MSA)

This example illustrates how to use Theseus as a Multiple Sequence Aligner. First, include the msa and general headers:
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
theseus::TheseusMSA aligner(penalties, heuristics, initial_sequence);
```

Once this is done, we can start adding sequences to our POA graph using the align functionality, that returns an Alignment object with CIGAR, path, and score information:
```
theseus::Alignment alignment_object = aligner.align(sequence);
```

Each time a new sequence is added to the POA graph, the graph is updated with the newly found variation (all the insertions, deletions and mismatches of the resulting alignment object).

Finally, we can output the result in four different formats: a graph in .gfa format, a multiple sequence alignment, a consensus sequence, and a graph in .dot format:
```
aligner.print_as_gfa(output_file);             // Output the compacted POA graph in .gfa format
aligner.print_as_msa(output_file);             // Output as a Multiple Sequence Alignment
sequence = aligner.get_consensus_sequence();   // Find the consensus sequence of the alignment
aligner.print_as_dot(output_file);             // Output the compacted POA graph in .dot format
```


### <a name="alignment_usage"></a> 2.2. Aligning a sequence to a graph

This example illustrates how to use Theseus as a general Sequence-to-Graph Aligner. First, include the alignment and general headers:
```
#include "theseus/alignment.h"
#include "theseus/graph.h"
#include "theseus/penalties.h"
#include "theseus/heuristics.h"
#include "theseus/theseus_aligner.h"
```

Then, create and configure an aligner object. This object is defined by three parameters: a set of penalties, a heuristics object, and a reference graph in Theseus' internal graph format. An example on how to set such parameters and create an aligner is found in the next code snippet:
```
theseus::Penalties penalties(match, mismatch, gap_open, gap_extend);
theseus::Heuristics heuristics(use_lag_pruning, use_density_drop);
theseus::TheseusAligner aligner(penalties, heuristics, std::move(graph));
```

Once this is done, we can start aligning sequences to the reference graph using the align functionality. **Importantly**, a call to the align function consists of four arguments: the **sequence** to be aligned, the **starting vertex** for the alignment, the **starting offset** in that starting vertex, and a boolean variable indicating the **direction of alignment**. The result of the alignment is an Alignment object with CIGAR, path and score information:
```
theseus::Alignment alignment_object = aligner.align(sequence, start_vertex, start_offset, align_reverse);
```


## <a name="theseus_tools"></a> 3.Tools

The Theseus library implements three minimal tools to use the Theseus algorithm on the MSA and Sequence-to-Graph alignment problems. It is important to note that these tools are not production ready.

### <a name="msa_tool"></a> 3.1. MSA tool: theseus_msa

This example illustrates how to use the **theseus_msa** tool. This tool computes the MSA of the set of sequences in an given input *.fasta* file. The executable is located in the path */build/tools/theseus_msa*:
```
cd build/tools/
```

Select the scoring scheme, set the input and output files and execute the tool. Each execution of *theseus_msa* lets you select the following parameters:
```
Usage: theseus_msa [OPTIONS]
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

An example of the execution of *theseus_msa* is shown in the following piece of code
```
./theseus_msa -m 0 -x 2 -o 3 -e 1 -t 0 -f output_file.out -s sequences.fasta
```


### <a name="seq_to_graph_tool"></a> 3.2. Seq-to-graph tool: theseus_aligner
This example illustrates how to use the **theseus_aligner** tool. This tool aligns a set of sequences, given their starting vertices and offsets, to a reference graph. Two input files are required: 1) The reference graph in *.gfa* format, 2) The sequences to be aligned and the starting alignment positions in *.fasta* format.

**[IMPORTANT]**
The *.fasta* file containing sequence and positional data has a special structure. As all .fasta files, the data associated to each sequence has two parts: 1) A line starting with ">" containing metadata, and 2) the sequence itself, that appears on the next lines.

1) Constists of three elements: *start_vertex* name, *start_offset* and orientation of the alignment. The orientation of the alignment can either be "+" or "-" and it indicates Theseus whether to align from the forward or reverse strand of the given start node. *[IMPORTANT]*, the *start_offset* value is respect the starting node in the given alignment orientation. That is, if alignment orientation is reverse ("-") and the offset is 0, alignment starts at the beggining of the reverse complemented segment.
2) Contains the sequence itself.

For instance, a sequence ACGT starting at vertex *v1*, offset *3* and in the forward orientation *+* would be codified as:
```
> v1 3 +
ACGT
> ... (following lines)
```

The executable is located in the path */build/tools/theseus_aligner*:
```
cd build/tools/
```

Select the scoring scheme, set the input and output files and execute the tool. Each execution of *theseus_aligner* lets you select the following parameters:
```
"Usage: theseus_aligner [OPTIONS]\n"
                 Options:\n"
                   Penalties:\n"
                   -m, --match <int>            The match penalty                                       [default=0]\n"
                   -x, --mismatch <int>         The mismatch penalty                                    [default=2]\n"
                   -o, --gapo <int>             The gap open penalty                                    [default=3]\n"
                   -e, --gape <int>             The gap extension penalty                               [default=1]\n\n"

                   I/O:\n"
                   -g, --graph_file <file>      Graph file in .gfa format                               [Required]\n"
                   -s, --sequences_file <file>  Sequences and starting positons in .fasta format        [Required]\n"
                   -f, --output_file <file>     Output file                                             [Required]\n\n"

                  Heuristics:\n"
                   -d  --density_heuristic     Activate the drop heuristic based on advancement density.            \n"
                   -l  --lag_pruning           Activate the pruning of diagonals lagging behind int the alignment.  \n";
```

An example of the execution of *theseus_aligner* is shown in the following piece of code
```
./theseus_aligner -m 0 -x 2 -o 3 -e 1 -g reference_graph.gfa -s sequences.fasta -f output.out
```

### <a name="consensus_tool"></a> 3.3. Consensus tool: pericles

This example illustrates how to use the **pericles** tool. This tool computes the MSA of the set of sequences in an given input *.fasta* file, allowing to add partial sequences, as long as they start on either end of a backbone sequence. The executable is located in the path */build/tools/pericles*:
```
cd build/tools/
```
**[IMPORTANT]**
The *.fasta* file containing sequences has a special structure. As all .fasta files, the data associated to each sequence has two parts: 1) A line starting with ">" containing metadata, and 2) the sequence itself, that appears on the next lines.

1) Constists of one element: ">" and a boolean value (0 or 1) indicating whether the associated read should be mapped in the canonical direction (from left to right and from source to sink) or in the reverse direction (from right to left and from sink to source).
2) Contains the sequence itself.

Examples aligning in the canonical and reversed directions:
```
>0
ACCT
>1
TCCA
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


## <a name="theseus_datasets"></a> 4. HEURISTICS
Theseus library implements some heuristic approaches that accelerate alignment at the expense of a limited loss in accuracy. In particular, Theseus implements 1) a **pruning heuristic** that discards diagonals that have fallen behind in the alignment, as long as the alignment has shown a significant advancement in the last scores, and 2) a **drop heuristic** that drops alignment when the advancement density (number of offsets advanced in the last scores) is very low.

You can activate these heuristics by initializing the Heuristics object with two separate boolean flags:
```
theseus::Heuristics heuristics(use_lag_pruning, use_density_drop);
```

Moreover, the three minimal tools provided in this alignment library allow you to activate these two heuristics from the command line, by adding the **-d** and **-l** flags:
```
./pericles -m 0 -x 2 -o 3 -e 1 -t 0 -f output_file.out -s sequences.fasta -d -l
```


## <a name="theseus_datasets"></a> 5. DATASETS

The datasets used in our paper are available in Zenodo. <!-- falta posar els links


## <a name="theseus_bugs"></a> . REPORTING BUGS AND FEATURE REQUEST

Feedback and bug reporting is highly appreciated. Please report any issue or suggestion on github or email to the main developer (albert.jimenez1@bsc.es). Don't hesitate to contact us if:
  - You experience any bug or crash.
  - You want to request a feature or have any suggestion.
  - Your application using the library is running slower than it should or you expected.
  - Need help integrating the library into your tool.


## <a name="theseus_licence"></a> 6. LICENSE

Theseus-lib is distributed under MIT licence.


## <a name="theseus_authors"></a> 7. AUTHORS

Albert Jimenez Blanco (albert.jimenez1@bsc.es) is the main developer and the person you should address your complaints.

Lorién López-Villellas has had major contributions both in the technical implementation of Theseus and the final structure of the library.

<!-- ## <a name="theseus_cite"></a> 7. CITATION

**Albert Jimenez-Blanco, Lorien Lopez-Villellas, Juan Carlos Moure, Miquel Moreto, Santiago Marco-Sola**. ["Theseus: Fast and Optimal Affine-Gap Sequence-to-Graph Alignment"](). Bioinformatics, 2026. -->

