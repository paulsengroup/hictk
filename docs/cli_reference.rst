..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

CLI Reference
#############

For an up-to-date list of subcommands and CLI options refer to hictk --help.

Subcommands
-----------

.. code-block:: text

  Blazing fast tools to work with .hic and .cool files.
  Usage: hictk [OPTIONS] SUBCOMMAND
  Options:
    -h,--help                   Print this help message and exit
    -V,--version                Display program version information and exit
  Subcommands:
    balance                     Balance Hi-C matrices using ICE, SCALE, or VC.
    convert                     Convert Hi-C matrices to a different format.
    dump                        Dump data from .hic and Cooler files to stdout.
    fix-mcool                   Fix corrupted .mcool files.
    load                        Build .cool and .hic files from interactions in various text formats.
    merge                       Merge multiple Cooler or .hic files into a single file.
    metadata                    Print file metadata to stdout.
    rename-chromosomes, rename-chroms
                                Rename chromosomes found in a Cooler file.
    validate                    Validate .hic and Cooler files.
    zoomify                     Convert single-resolution Cooler and .hic files to multi-resolution by coarsening.

hictk balance
-------------

.. code-block:: text

  Balance Hi-C matrices using ICE, SCALE, or VC.
  Usage: hictk balance [OPTIONS] [SUBCOMMAND]
  Options:
    -h,--help                   Print this help message and exit
  Subcommands:
    ice                         Balance Hi-C matrices using ICE.
    scale                       Balance Hi-C matrices using SCALE.
    vc                          Balance Hi-C matrices using VC.

hictk balance ice
-----------------

.. code-block:: text

  Balance Hi-C matrices using ICE.
  Usage: hictk balance ice [OPTIONS] input
  Positionals:
    input TEXT:((HiC) OR (Cooler)) OR (Multires-cooler) REQUIRED
                                Path to the .hic, .cool or .mcool file to be balanced.
  Options:
    -h,--help                   Print this help message and exit
    --mode TEXT:{gw,trans,cis} [gw]
                                Balance matrix using:
                                 - genome-wide interactions (gw)
                                 - trans-only interactions (trans)
                                 - cis-only interactions (cis)
    --tmpdir TEXT:DIR           Path to a folder where to store temporary data.
    --ignore-diags UINT [2]     Number of diagonals (including the main diagonal) to mask before balancing.
    --mad-max FLOAT:NONNEGATIVE [5]
                                Mask bins using the MAD-max filter.
                                bins whose log marginal sum is less than --mad-max median
                                absolute deviations below the median log marginal sum of
                                all the bins in the same chromosome.
    --min-nnz UINT [10]         Mask rows with fewer than --min-nnz non-zero entries.
    --min-count UINT [0]        Mask rows with fewer than --min-count interactions.
    --tolerance FLOAT:NONNEGATIVE [1e-05]
                                Threshold of the variance of marginals used to determine whether
                                the algorithm has converged.
    --max-iters UINT:POSITIVE [500]
                                Maximum number of iterations.
    --rescale-weights,--no-rescale-weights{false}
                                Rescale weights such that rows sum approximately to 2.
    --name TEXT                 Name to use when writing weights to file.
                                Defaults to ICE, INTER_ICE and GW_ICE when --mode is cis, trans and gw, respectively.
    --create-weight-link,--no-create-weight-link{false}
                                Create a symbolic link to the balancing weights at clr::/bins/weight.
                                Ignored when balancing .hic files
    --in-memory                 Store all interactions in memory (greatly improves performance).
    --stdout                    Write balancing weights to stdout instead of writing them to the input file.
    --chunk-size UINT:POSITIVE [10000000]
                                Number of interactions to process at once. Ignored when using --in-memory.
    -v,--verbosity UINT:INT in [1 - 4] []
                                Set verbosity of output to the console.
    -t,--threads UINT:UINT in [1 - 32] [1]
                                Maximum number of parallel threads to spawn.
    -l,--compression-lvl UINT:INT in [0 - 19] []
                                Compression level used to compress temporary files using ZSTD.
    -f,--force                  Overwrite existing files and datasets (if any).

hictk balance scale
-------------------

.. code-block:: text

  Balance Hi-C matrices using SCALE.
  Usage: hictk balance scale [OPTIONS] input
  Positionals:
    input TEXT:((HiC) OR (Cooler)) OR (Multires-cooler) REQUIRED
                                Path to the .hic, .cool or .mcool file to be balanced.
  Options:
    -h,--help                   Print this help message and exit
    --mode TEXT:{gw,trans,cis} [gw]
                                Balance matrix using:
                                 - genome-wide interactions (gw)
                                 - trans-only interactions (trans)
                                 - cis-only interactions (cis)
    --tmpdir TEXT               Path to a folder where to store temporary data.
    --max-percentile FLOAT [10]
                                Percentile used to compute the maximum number of nnz values that cause a row to be masked.
    --max-row-sum-err FLOAT:NONNEGATIVE [0.05]
                                Row sum threshold used to determine whether convergence has been achieved.
    --tolerance FLOAT:NONNEGATIVE [0.0001]
                                Threshold of the variance of marginals used to determine whether
                                the algorithm has converged.
    --max-iters UINT:POSITIVE [500]
                                Maximum number of iterations.
    --rescale-weights,--no-rescale-weights{false}
                                Rescale weights such that the sum of the balanced matrix is similar
                                to that of the input matrix.
    --name TEXT                 Name to use when writing weights to file.
                                Defaults to SCALE, INTER_SCALE and GW_SCALE when --mode is cis, trans and gw, respectively.
    --create-weight-link,--no-create-weight-link{false}
                                Create a symbolic link to the balancing weights at clr::/bins/weight.
                                Ignored when balancing .hic files
    --in-memory                 Store all interactions in memory (greatly improves performance).
    --stdout                    Write balancing weights to stdout instead of writing them to the input file.
    --chunk-size UINT:POSITIVE [10000000]
                                Number of interactions to process at once. Ignored when using --in-memory.
    -v,--verbosity UINT:INT in [1 - 4] []
                                Set verbosity of output to the console.
    -t,--threads UINT:UINT in [1 - 32] [1]
                                Maximum number of parallel threads to spawn.
    -l,--compression-lvl UINT:INT in [0 - 19] []
                                Compression level used to compress temporary files using ZSTD.
    -f,--force                  Overwrite existing files and datasets (if any).

hictk balance vc
----------------

.. code-block:: text

  Balance Hi-C matrices using VC.
  Usage: hictk balance vc [OPTIONS] input
  Positionals:
    input TEXT:((HiC) OR (Cooler)) OR (Multires-cooler) REQUIRED
                                Path to the .hic, .cool or .mcool file to be balanced.
  Options:
    -h,--help                   Print this help message and exit
    --mode TEXT:{gw,trans,cis} [gw]
                                Balance matrix using:
                                 - genome-wide interactions (gw)
                                 - trans-only interactions (trans)
                                 - cis-only interactions (cis)
    --rescale-weights,--no-rescale-weights{false}
                                Rescale weights such that the sum of the balanced matrix is similar
                                to that of the input matrix.
    --name TEXT                 Name to use when writing weights to file.
                                Defaults to VC, INTER_VC and GW_VC when --mode is cis, trans and gw, respectively.
    --create-weight-link,--no-create-weight-link{false}
                                Create a symbolic link to the balancing weights at clr::/bins/weight.
                                Ignored when balancing .hic files
    --stdout                    Write balancing weights to stdout instead of writing them to the input file.
    -v,--verbosity UINT:INT in [1 - 4] []
                                Set verbosity of output to the console.
    -f,--force                  Overwrite existing files and datasets (if any).

hictk convert
-------------

.. code-block:: text

  Convert Hi-C matrices to a different format.
  Usage: hictk convert [OPTIONS] input output
  Positionals:
    input TEXT:((HiC) OR (Cooler)) OR (Multires-cooler) REQUIRED
                                Path to the .hic, .cool or .mcool file to be converted.
    output TEXT REQUIRED        Output path. File extension is used to infer output format.
  Options:
    -h,--help                   Print this help message and exit
    --output-fmt TEXT:{cool,mcool,hic} [auto]
                                Output format (by default this is inferred from the output file extension).
                                Should be one of:
                                - cool
                                - mcool
                                - hic
    -r,--resolutions UINT:POSITIVE ...
                                One or more resolutions to be converted. By default all resolutions are converted.
    --normalization-methods TEXT [ALL]  ...
                                Name of one or more normalization methods to be copied.
                                By default, vectors for all known normalization methods are copied.
                                Pass NONE to avoid copying normalization vectors.
    --fail-if-norm-not-found    Fail if any of the requested normalization vectors are missing.
    -g,--genome TEXT            Genome assembly name. By default this is copied from the .hic file metadata.
    --tmpdir TEXT:DIR           Path where to store temporary files.
    --chunk-size UINT:POSITIVE [10000000]
                                Batch size to use when converting .[m]cool to .hic.
    -v,--verbosity UINT:INT in [1 - 4] []
                                Set verbosity of output to the console.
    -t,--threads UINT:UINT in [2 - 12] [2]
                                Maximum number of parallel threads to spawn.
                                When converting from hic to cool, only two threads will be used.
    -l,--compression-lvl UINT:INT in [1 - 32] [6]
                                Compression level used to compress interactions.
                                Defaults to 6 and 10 for .cool and .hic files, respectively.
    --skip-all-vs-all,--no-skip-all-vs-all{false}
                                Do not generate All vs All matrix.
                                Has no effect when creating .[m]cool files.
    --count-type TEXT:{auto,int,float} [auto]
                                Specify the strategy used to infer count types when converting .hic
                                files to .[m]cool format.
                                Can be one of: int, float, auto.
    -f,--force                  Overwrite existing files (if any).

hictk dump
----------

.. code-block:: text

  Dump data from .hic and Cooler files to stdout.
  Usage: hictk dump [OPTIONS] uri
  Positionals:
    uri TEXT:(((HiC) OR (Cooler)) OR (Multires-cooler)) OR (Single-cell-cooler) REQUIRED
                                Path to a .hic, .cool or .mcool file (Cooler URI syntax supported).
  Options:
    -h,--help                   Print this help message and exit
    --resolution UINT:NONNEGATIVE
                                Hi-C matrix resolution (ignored when file is in .cool format).
    --matrix-type ENUM:value in {expected->2,observed->0,oe->1} OR {2,0,1} [observed]
                                Matrix type (ignored when file is not in .hic format).
    --matrix-unit ENUM:value in {BP->0,FRAG->1} OR {0,1} [BP]
                                Matrix unit (ignored when file is not in .hic format).
    -t,--table TEXT:{chroms,bins,pixels,normalizations,resolutions,cells,weights} [pixels]
                                Name of the table to dump.
    -r,--range TEXT [all]  Excludes: --query-file --cis-only --trans-only
                                Coordinates of the genomic regions to be dumped following UCSC-style notation (chr1:0-1000).
    --range2 TEXT [all]  Needs: --range Excludes: --query-file --cis-only --trans-only
                                Coordinates of the genomic regions to be dumped following UCSC-style notation (chr1:0-1000).
    --query-file TEXT:(FILE) OR ({-}) Excludes: --range --range2 --cis-only --trans-only
                                Path to a BEDPE file with the list of coordinates to be fetched (pass - to read queries from stdin).
    --cis-only Excludes: --range --range2 --query-file --trans-only
                                Dump intra-chromosomal interactions only.
    --trans-only Excludes: --range --range2 --query-file --cis-only
                                Dump inter-chromosomal interactions only.
    -b,--balance TEXT [NONE]    Balance interactions using the given method.
    --sorted,--unsorted{false}  Return interactions in ascending order.
    --join,--no-join{false}     Output pixels in BG2 format.

hictk fix-mcool
---------------

.. code-block:: text

  Fix corrupted .mcool files.
  Usage: hictk fix-mcool [OPTIONS] input output
  Positionals:
    input TEXT:Multires-cooler REQUIRED
                                Path to a corrupted .mcool file.
    output TEXT REQUIRED        Path where to store the restored .mcool.
  Options:
    -h,--help                   Print this help message and exit
    --tmpdir TEXT:DIR           Path to a folder where to store temporary data.
    --skip-balancing            Do not recompute or copy balancing weights.
    --check-base-resolution     Check whether the base resolution is corrupted.
    --in-memory                 Store all interactions in memory while balancing (greatly improves performance).
    --chunk-size UINT:POSITIVE [10000000]
                                Number of interactions to process at once during balancing.
                                Ignored when using --in-memory.
    -v,--verbosity UINT:INT in [1 - 4] []
                                Set verbosity of output to the console.
    -t,--threads UINT:UINT in [1 - 32] [1]
                                Maximum number of parallel threads to spawn (only applies to the balancing stage).
    -l,--compression-lvl UINT:INT in [0 - 19] []
                                Compression level used to compress temporary files using ZSTD (only applies to the balancing stage).
    -f,--force                  Overwrite existing files (if any).

hictk load
----------

.. code-block:: text

  Build .cool and .hic files from interactions in various text formats.
  Usage: hictk load [OPTIONS] interactions output-path
  Positionals:
    interactions TEXT:(FILE) OR ({-}) REQUIRED
                                Path to a file with the interactions to be loaded.
                                Common compression formats are supported (namely, gzip, bzip2, lz4, lzo, and zstd).
                                Pass "-" to indicate that interactions should be read from stdin.
    output-path TEXT REQUIRED   Path to output file.
                                File extension will be used to infer the output format.
                                This behavior can be overridden by explicitly specifying an
                                output format through option --output-fmt.
  Options:
    -h,--help                   Print this help message and exit
    -c,--chrom-sizes TEXT:FILE Excludes: --bin-table
                                Path to .chrom.sizes file.
                                Required when interactions are not in 4DN pairs format.
    -b,--bin-size UINT:POSITIVE Excludes: --bin-table
                                Bin size (bp).
                                Required when --bin-table is not used.
    --bin-table TEXT:FILE Excludes: --chrom-sizes --bin-size
                                Path to a BED3+ file with the bin table.
    -f,--format TEXT:{4dn,validpairs,bg2,coo} REQUIRED
                                Input format.
    --output-fmt TEXT:{auto,cool,hic} [auto]
                                Output format (by default this is inferred from the output file extension).
                                Should be one of:
                                - autol
                                - cool
                                - hic
    --force                     Force overwrite existing output file(s).
    --assembly TEXT [unknown]   Assembly name.
    --drop-unknown-chroms       Ignore records referencing unknown chromosomes.
    --one-based,--zero-based{false}
                                Interpret genomic coordinates or bins as one/zero based.
                                By default coordinates are assumed to be one-based for interactions in
                                4dn and validpairs formats and zero-based otherwise.
    --count-as-float            Interactions are floats.
    --skip-all-vs-all,--no-skip-all-vs-all{false}
                                Do not generate All vs All matrix.
                                Has no effect when creating .cool files.
    --assume-sorted,--assume-unsorted{false}
                                Assume input files are already sorted.
    --chunk-size UINT [10000000]
                                Number of pixels to buffer in memory.
    -l,--compression-lvl UINT:INT bounded to [1 - 32]
                                Compression level used to compress interactions.
                                Defaults to 6 and 10 for .cool and .hic files, respectively.
    -t,--threads UINT:UINT in [2 - 12] [2]
                                Maximum number of parallel threads to spawn.
                                When loading interactions in a .cool file, only up to two threads will be used.
    --tmpdir TEXT:DIR           Path to a folder where to store temporary data.
    -v,--verbosity UINT:INT in [1 - 4] []
                                Set verbosity of output to the console.

hictk merge
-----------

.. code-block:: text

  Merge multiple Cooler or .hic files into a single file.
  Usage: hictk merge [OPTIONS] input-files...
  Positionals:
    input-files TEXT:(Cooler) OR (HiC) x 2 REQUIRED
                                Path to two or more Cooler or .hic files to be merged (Cooler URI syntax supported).
  Options:
    -h,--help                   Print this help message and exit
    -o,--output-file TEXT REQUIRED
                                Output Cooler or .hic file (Cooler URI syntax supported).
    --output-fmt TEXT:{cool,hic} [auto]
                                Output format (by default this is inferred from the output file extension).
                                Should be one of:
                                - cool
                                - hic
    --resolution UINT:NONNEGATIVE
                                Hi-C matrix resolution (ignored when input files are in .cool format).
    -f,--force                  Force overwrite output file.
    --chunk-size UINT [10000000]
                                Number of pixels to store in memory before writing to disk.
    -l,--compression-lvl UINT:INT bounded to [1 - 32]
                                Compression level used to compress interactions.
                                Defaults to 6 and 10 for .cool and .hic files, respectively.
    -t,--threads UINT:UINT in [1 - 32] [1]
                                Maximum number of parallel threads to spawn.
                                When merging interactions in Cooler format, only a single thread will be used.
    --tmpdir TEXT:DIR           Path to a folder where to store temporary data.
    --skip-all-vs-all,--no-skip-all-vs-all{false}
                                Do not generate All vs All matrix.
                                Has no effect when merging .cool files.
    --count-type TEXT:{int,float} [int]
                                Specify the count type to be used when merging files.
                                Ignored when the output file is in .hic format.
    -v,--verbosity UINT:INT in [1 - 4] []
                                Set verbosity of output to the console.

hictk metadata
--------------

.. code-block:: text

  Print file metadata to stdout.
  Usage: hictk metadata [OPTIONS] uri
  Positionals:
    uri TEXT:(((Cooler) OR (Multires-cooler)) OR (Single-cell-cooler)) OR (HiC) REQUIRED
                                Path to a .hic or .[ms]cool file (Cooler URI syntax supported).
  Options:
    -h,--help                   Print this help message and exit
    -f,--output-format TEXT:{json,toml,yaml} [json]
                                Format used to return file metadata.
                                Should be one of: json, toml, or yaml.
    --include-file-path,--exclude-file-path{false}
                                Output the given input path using attribute "uri".
    --recursive                 Print metadata for each resolution or cell contained in a
                                multi-resolution or single-cell file.

hictk rename-chromosomes
------------------------

.. code-block:: text

  Rename chromosomes found in a Cooler file.
  Usage: hictk rename-chromosomes [OPTIONS] uri
  Positionals:
    uri TEXT REQUIRED           Path to a or .[ms]cool file (Cooler URI syntax supported).
  Options:
    -h,--help                   Print this help message and exit
    --name-mappings TEXT Excludes: --add-chr-prefix --remove-chr-prefix
                                Path to a two column TSV with pairs of chromosomes to be renamed.
                                The first column should contain the original chromosome name,
                                while the second column should contain the destination name to use when renaming.
    --add-chr-prefix Excludes: --name-mappings --remove-chr-prefix
                                Prefix chromosome names with "chr".
    --remove-chr-prefix Excludes: --name-mappings --add-chr-prefix
                                Remove prefix "chr" from chromosome names.
    -v,--verbosity UINT:INT in [1 - 4] []
                                Set verbosity of output to the console.

hictk validate
--------------

.. code-block:: text

  Validate .hic and Cooler files.
  Usage: hictk validate [OPTIONS] uri
  Positionals:
    uri TEXT REQUIRED           Path to a .hic or .[ms]cool file (Cooler URI syntax supported).
  Options:
    -h,--help                   Print this help message and exit
    --validate-index            Validate Cooler index (may take a long time).
    -f,--output-format TEXT:{json,toml,yaml} [json]
                                Format used to report the outcome of file validation.
                                Should be one of: json, toml, or yaml.
    --include-file-path,--exclude-file-path{false}
                                Output the given input path using attribute "uri".
    --exhaustive,--fail-fast{false}
                                When processing multi-resolution or single-cell files,
                                do not fail as soon as the first error is detected.
    --quiet                     Don't print anything to stdout. Success/failure is reported through exit codes

hictk zoomify
-------------

.. code-block:: text

  Convert single-resolution Cooler and .hic files to multi-resolution by coarsening.
  Usage: hictk zoomify [OPTIONS] cooler/hic [m]cool/hic
  Positionals:
    cooler/hic TEXT:(Cooler) OR (HiC) REQUIRED
                                Path to a .cool or .hic file (Cooler URI syntax supported).
    [m]cool/hic TEXT REQUIRED   Output path.
                                When zoomifying Cooler files, providing a single resolution through
                                --resolutions and specifying --no-copy-base-resolution, the output file
                                will be in .cool format.
  Options:
    -h,--help                   Print this help message and exit
    --force                     Force overwrite existing output file(s).
    --resolutions UINT ...      One or more resolutions to be used for coarsening.
    --copy-base-resolution,--no-copy-base-resolution{false}
                                Copy the base resolution to the output file.
    --nice-steps,--pow2-steps{false} [--nice-steps]
                                Use nice or power of two steps to automatically generate the list of resolutions.
                                Example:
                                Base resolution: 1000
                                Pow2: 1000, 2000, 4000, 8000...
                                Nice: 1000, 2000, 5000, 10000...
    -l,--compression-lvl UINT:INT bounded to [1 - 32] [6]
                                Compression level used to compress interactions.
                                Defaults to 6 and 10 for .mcool and .hic files, respectively.
    -t,--threads UINT:UINT in [1 - 32] [1]
                                Maximum number of parallel threads to spawn.
                                When zoomifying interactions from a .cool file, only a single thread will be used.
    --chunk-size UINT [10000000]
                                Number of pixels to buffer in memory.
                                Only used when zoomifying .hic files.
    --skip-all-vs-all,--no-skip-all-vs-all{false}
                                Do not generate All vs All matrix.
                                Has no effect when zoomifying .cool files.
    --tmpdir TEXT:DIR           Path to a folder where to store temporary data.
    -v,--verbosity UINT:INT in [1 - 4] []
                                Set verbosity of output to the console.
