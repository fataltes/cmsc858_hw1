Succinct Data Structures and Operations
========================================

# Introduction
Succinct operations, specifically rank and select over bit-vectors are important in helping with time and space efficiency 
of memory and time critical softwares working with terabytes of data.
They operate over a set of bits.For a bit-vector *bv* with size *bv.size*, *Rank1(idx)* for *0<=idx<bv.size*
returns the number of *1*s until index `idx` inclusive and *Rank0(idx)* returns number of *0*s in the same range.
Select operation can be counted as the inverse of rank where *select1(v)* for *1<=v<=Rank1(bv.size)* returns the index of the *v*th *1* in *bv*.
To apply rank and select operations in a world with alphabets more than only *0* and *1*, one of the very popular
and practical data structures are the Wavelet trees. Wavelet trees are binary perfectly balanced trees
that index a text consisting of alphabet in set $\sigma$ to a set of *0* and *1*s so that operations such as rank of the character *c*
and select of character *c*, $c \in \sigma$, can be reduced to a bunch of bit-vector rank and select operations.


In this homework, we were supposed to implement rank and select operations for an already existing bit-vector in tasks 1 and 2
and later, use these operations in querying a wavelet tree data structure that we have constructed, stored to, and later loaded from the disk.

All code has been implemented in C++14 and is available in the CMake-based project repository "*https://github.com/fataltes/cmsc858_hw1.git*".
I have used a slightly modified version of the `compact bit-vector` package from "*https://github.com/gmarcais/compact_vector.git*" as the underlying bit-vector data structure.

This report only contains the abstract overview of the implementation and the reports and plots.

The installation process is the typical cmake, build process. The least required cmake version is *3.9*
and I have compiled all the codes using gcc version *5.5*.

The main command available is `>bvOperators` and later we can add options and sub-commands to see the results for different tasks.

Note: In general, the main challenge in coding in C++ is handling memory, being careful about bounds, dangling pointers, bad allocs, segfaults, and things like that.

I think initializing a project from scratch in C++ which should be a reasonably well-structures and working one was also another interesting challenge.


## Task1
In this task, I implemented a rank data structure over the compact::bitvector structure.
The rank data structure itself consists of three bit-vectors, R_s, R_b, and R_p. I tried the version that
trades off time for space? so I will store R_p as well which is $\bigO(log(log(n))*log(n)*2^log(n))$ extra space,
but reduces the time to just the constant time required to have a couple look-ups over the three bit-vectors
rather than at the point of the blocks, walking them until we find the desired value/index.
The following are the average time required per rank query for various sizes of bit-vectors,
where *1/10*th of each bit-vector bits are set and also the overhead size of the rank data structure for each sample.


plots go here



The most difficult section in this task was undoubtedly handling R_p indices and later finding the appropriate index
for each bit-vector of R_s, R_b, and R_p based on the input requested value.

You can run the command for this report as following:
```
> ./bvOperators report -t rank -p console
```
If you want to store the results in a file, provide the parent directory path after `-p` rather than the keyword `console'.
```
> ./bvOperators report -t rank -p <par_dir>
```
There are 3 options available for report sub-command to set minimum, maximum, and jump size for different bit-vectors that rank is going to be tested on.

# Task2
In this task, I implemented the *O(log(n))* time version of the select which requires no extra space
if rank is also available on the bit-vector and the same space as rank otherwise. It basically is
a binary search using rank operation to find the bit with the appropriate rank and return its index.

There was not any main challenges in this task except maybe in the engineering aspect of it.
The most challenging part for me was how to prevent code chunk repetition by best utilizing templates
etc.

Same plots as task one have been provided for this task as well.


plots go here



You can run the command for this report as following:
```
> ./bvOperators report -t select -p console
```
If you want to store the results in a file, provide the parent directory path after `-p` rather than the keyword `console`.
```
> ./bvOperators report -t select -p <par_dir>
```
There are 3 options available for report sub-command to set minimum, maximum, and jump size for different bit-vectors that rank is going to be tested on.

# Task3
In this section, I constructed the wavelet tree based on the serialized algorithm that was explained in the paper
``Simple, Fast and Lightweight Parallel Wavelet Tree Construction'' and stored it in the Level-Wise Wavelet Tree data structure
which is a bit-vector containing all the bits required to store different levels of the tree.
What I store at the end is the wavelet tree bit-vector in binary format, count of unique characters, list of the characters in sorted order,
each of a size of byte and the original sequence length in number of characters.

The construction Memory for my algorithm small and limited to the fixed buffer size that I read from the file and put there
plus the start position vector. I don't load the whole file all at once into memory and read and process it buffer by buffer.
For that reason, I have two passes over the file, but based on a few engineering tweaks I managed to only read the file twice and
construct the whole tree. I also don't store the last level of the tree which is all *0*s as each block
belongs to one single character (or no characters) at that level and hence there is no distinction required between
characters ending up in left and right branch. I also don't store the block start positions in file as it can be
reconstructed easily every time that we load the index.

For improving query time (access, rank, and select) every time that I load the index in addition to constructing
the vector of start positions for each block, I also keep a vector of rank of *1* at the start position of
each block. That saves a bunch of bitwise rank queries required in each of the character query processes
along the tree and although it does not affect the asymptotics of the query time, it helps reduce the constant
value and make the processes faster in trade for a small amount of memory. There is one more interesting 
engineering I did for improving the select performance on the wavelet tree too.
At each level, when we want to call a select which later calls a binary search using rank, 
I instead directly call the binary search with the manual start and end positions
 being the start and end position for the block I'm looking into.
 In this way, at each level, rather than doing a binary search over a sequence of size n
 which in total costs $O(nlog(\sigma))$,
 I do the search over a sequence of approximate length ~$n/2^level$. 

The trickiest part in this task for me was understanding and later implementing the process of
selecting the smallest index `idx` that a character *c* occurs *s* time up to that index.
It was challenging because unlike rank that narrows as we go down the levels of the tree
and makes it possible to walk the tree top-down, in select, we need to start from the lowest level
which is actually the level that each character belongs to a single block and then walk up the tree.

Here are the plots showing the scalability of the data structure with respect to the bit-vector size
and ....


Since I am putting all of these as one executable, the commands interfaces are a little bit
different from what the professor has specifically asked but they answer the same queries and
build process that is required in the homework.
In the following lines, I bring a few examples of different sub-commands for constructing and querying the wavelet tree.

To construct the wavelet tree, we require two inputs, the input file containing the sequence `<input>`
and the address pointing to the directory we will store the index in, `<idx_par_dir>`.

```
> ./bvOperators wv_build -i <input> -p <idx_par_dir>
```

To call different queries of <access, rank, and select> on the wavelet tree,
we require two main inputs, the query file `<qfile>` and the index directory `<idx_par_dir>`
in addition to specifying the type of the query.[^ps]
```
> ./bvOperators wv -t <access/rank/select> -p <par_dir> -i <qfile>
```



---------------
[^ps]: P.S. Although I tried to catch a couple special cases with an invalid input
such as querying for a character that does not exist in the alphabet of wavelet tree,
Unfortunately, for the sake of time-limitation for this homework I was not able to manage all the
specific flaws in the input arguments such as a wrong file format, a non-existing directory, etc.
However, these verifications are in the todo list of this project in the future.