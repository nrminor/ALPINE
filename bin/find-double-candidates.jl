#!/usr/bin/env julia

# load pipeline module
using LongInfectionFinder

# load I/O
const high_distance_metadata = islink(ARGS[1]) ? readlink(ARGS[1]) : ARGS[1]
const anachronistic_metadata = islink(ARGS[2]) ? readlink(ARGS[2]) : ARGS[2]
const anachronistic_sequences = islink(ARGS[3]) ? readlink(ARGS[3]) : ARGS[3]

find_double_candidates(high_distance_metadata, anachronistic_metadata, anachronistic_sequences)
