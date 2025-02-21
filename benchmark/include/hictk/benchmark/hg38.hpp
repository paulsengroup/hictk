// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <vector>

#include "hictk/chromosome.hpp"

namespace hictk::benchmark {

// clang-format off
// NOLINTNEXTLINE(cert-err58-cpp)
inline const std::vector hg38{
        hictk::Chromosome{0,  "chr1",  248956422},
        hictk::Chromosome{1,  "chr2",  242193529},
        hictk::Chromosome{2,  "chr3",  198295559},
        hictk::Chromosome{3,  "chr4",  190214555},
        hictk::Chromosome{4,  "chr5",  181538259},
        hictk::Chromosome{5,  "chr6",  170805979},
        hictk::Chromosome{6,  "chr7",  159345973},
        hictk::Chromosome{7,  "chr8",  145138636},
        hictk::Chromosome{8,  "chr9",  138394717},
        hictk::Chromosome{9,  "chr10", 133797422},
        hictk::Chromosome{10, "chr11", 135086622},
        hictk::Chromosome{11, "chr12", 133275309},
        hictk::Chromosome{12, "chr13", 114364328},
        hictk::Chromosome{13, "chr14", 107043718},
        hictk::Chromosome{14, "chr15", 101991189},
        hictk::Chromosome{15, "chr16", 90338345},
        hictk::Chromosome{16, "chr17", 83257441},
        hictk::Chromosome{17, "chr18", 80373285},
        hictk::Chromosome{18, "chr19", 58617616},
        hictk::Chromosome{19, "chr20", 64444167},
        hictk::Chromosome{20, "chr21", 46709983},
        hictk::Chromosome{21, "chr22", 50818468},
        hictk::Chromosome{22, "chrX",  156040895},
        hictk::Chromosome{23, "chrY",  57227415}
};
// clang-format on

}  // namespace hictk::benchmark
