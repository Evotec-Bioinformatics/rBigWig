x <- rBigWig::fetch_region("/data/projects/targetbcd/Alignments/SE2002_1_CSFP200001243-1a_HT7GYDSXX_L1_1_bwa/SE2002_1_CSFP200001243-1a_HT7GYDSXX_L1_1.bw",
                           "21", 1000, 10000000)


x <- rBigWig::fetch_region_means("/data/projects/targetbcd/Alignments/SE2002_1_CSFP200001243-1a_HT7GYDSXX_L1_1_bwa/SE2002_1_CSFP200001243-1a_HT7GYDSXX_L1_1.bw",
                           "21", 1000, 10000000, 100)
