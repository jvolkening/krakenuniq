import re

def edit(path, pairs, count=1):
    s = open(path).read()
    for a, b in pairs:
        assert s.count(a) == count, (path, a, s.count(a))
        s = s.replace(a, b)
    open(path, 'w').write(s)

# krakendb.hpp
edit('src/krakendb.hpp', [
    ("void prepare_chunking(const uint64_t max_bytes_for_db);",
     "// If shared_buffer is given (at least max_bytes_for_db bytes, owned by the caller),\n"
     "    // it is used instead of allocating a private buffer.\n"
     "    void prepare_chunking(const uint64_t max_bytes_for_db, char *shared_buffer = NULL);"),
    ("    char *data;\n    uint64_t data_size;\n    uint64_t data_offset;\n  };\n}",
     "    char *data;\n    bool owns_data;\n    uint64_t data_size;\n    uint64_t data_offset;\n  };\n}"),
])

# krakendb.cpp
edit('src/krakendb.cpp', [
    ("  data = NULL;\n  data_size = 0;\n", "  data = NULL;\n  owns_data = false;\n  data_size = 0;\n"),
], count=2)
s = open('src/krakendb.cpp').read()
a = s.index("void KrakenDB::prepare_chunking")
b = s.index("  this->_chunks = 0;\n\n  const uint64_t max_idx_pos")
new = '''void KrakenDB::prepare_chunking(const uint64_t max_bytes_for_db, char *shared_buffer) {
  if (shared_buffer != NULL) {
    data = shared_buffer;
    owns_data = false;
  } else {
    data = (char*) mmap(0, max_bytes_for_db, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (data == (void *) -1)
    {
      printf("ERROR: Could not allocated memory for loading database chunks (error code %i: %s)\\n", errno, strerror(errno));
      exit(-1);
    }
    owns_data = true;
    mlock(data, max_bytes_for_db);
  }
  data_size = max_bytes_for_db;

'''
s = s[:a] + new + s[b:]
open('src/krakendb.cpp', 'w').write(s)
edit('src/krakendb.cpp', [
    ("if (data != NULL && data_size > 0)\n    munmap", "if (owns_data && data != NULL && data_size > 0)\n    munmap"),
    ("    idx_pos = upper_bound(idx_pos, max_idx_pos);\n",
     "    const uint64_t prev_idx_pos = idx_pos;\n"
     "    idx_pos = upper_bound(idx_pos, max_idx_pos);\n"
     "    if (idx_pos == prev_idx_pos) {\n"
     "      fprintf(stderr, \"ERROR: The preload size (%llu bytes) is too small to hold a single minimizer bin of the database.\\n\", (unsigned long long) max_bytes_for_db);\n"
     "      exit(1);\n"
     "    }\n"),
])

# classify.cpp
edit('src/classify.cpp', [
    ("void merge_intermediate_results_by_workers(const bool first_intermediate_output, const std::string & tmp_file_name) {",
     "// same_db_as_prev: the previous intermediate result came from another chunk of the same database\n"
     "void merge_intermediate_results_by_workers(const bool first_intermediate_output, const bool same_db_as_prev, const std::string & tmp_file_name) {"),
    ("        assert(!(taxa[i] && taxa_prev_summary[i])); // a kmer cannot be found in more than one database chunk",
     "        // a kmer cannot be found in more than one chunk of the same database (but may be in several databases)\n"
     "        assert(!same_db_as_prev || !(taxa[i] && taxa_prev_summary[i]));"),
    ("          taxa[i] = taxa_prev_summary[i]; // if previous chunk/database got a match, keep the tax info",
     "          taxa[i] = taxa_prev_summary[i]; // if previous chunk/database got a match, keep the tax info (earlier databases take precedence)"),
    ("merge_intermediate_results_by_workers(first_intermediate_output, tmp_file_name);",
     "merge_intermediate_results_by_workers(first_intermediate_output, db_chunk_id > 0, tmp_file_name);"),
    ("db_chunk_id < KrakenDatabases[0]->chunks();", "db_chunk_id < KrakenDatabases[i]->chunks();"),
    ("  static vector<KrakenDBIndex> db_indices (DB_filenames.size());\n",
     "  static vector<KrakenDBIndex> db_indices (DB_filenames.size());\n\n"
     "  // A single buffer shared by all databases (they are processed one after the other),\n"
     "  // so that the memory used for chunks never exceeds Populate_memory_size\n"
     "  char *chunk_buffer = NULL;\n"
     "  if (Populate_memory && Populate_memory_size > 0)\n"
     "  {\n"
     "    chunk_buffer = (char*) mmap(0, Populate_memory_size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);\n"
     "    if (chunk_buffer == (void *) -1)\n"
     "    {\n"
     "      cerr << \"ERROR: Could not allocate \" << Populate_memory_size << \" bytes for loading database chunks: \" << strerror(errno) << endl;\n"
     "      return 1;\n"
     "    }\n"
     "    mlock(chunk_buffer, Populate_memory_size);\n"
     "  }\n"),
    ("prepare_chunking(Populate_memory_size);", "prepare_chunking(Populate_memory_size, chunk_buffer);"),
])

# scripts/krakenuniq
s = open('scripts/krakenuniq').read()
a = s.index("if (scalar(@db_prefix) > 1 && defined($preload_size) ) {")
b = s.index("}\n", s.index("exit 1;", a)) + 2
s = s[:a] + s[b:].lstrip("\n")
old = " overrides --preload flag.  Can only be used with a single input database."
assert s.count(old) == 1
s = s.replace(old, " overrides --preload flag.  With multiple databases, SIZE is the total memory used for all of them.")
open('scripts/krakenuniq', 'w').write(s)

# README.md
edit('README.md', [
    ("IMPORTANT!  The --preload-size option can only be used with a single input database.",
     "With several input databases (hierarchical classification), the chunks of all databases are loaded one after the other into the same buffer, so the amount given to --preload-size is the total memory used for the databases. The results are identical to those without --preload-size. The amount must be large enough to hold the largest single chunk (one minimizer bin plus its index entries) of any database; otherwise an error is reported."),
    ("overrides --preload flag\n", "overrides --preload flag. With multiple databases, SIZE is the total memory used for all of them\n"),
])
