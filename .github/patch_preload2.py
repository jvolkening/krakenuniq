def edit(path, pairs):
    s = open(path).read()
    for a, b in pairs:
        assert s.count(a) == 1, (path, a, s.count(a))
        s = s.replace(a, b)
    open(path, 'w').write(s)

# query the database that the current chunk belongs to (was hard-coded to database 0)
edit('src/classify.cpp', [
    ("classify_sequence_with_db_chunk(work_unit[j], fp, db_chunk_id, 0);",
     "classify_sequence_with_db_chunk(work_unit[j], fp, db_chunk_id, i);"),
])

# bins in the tiny test databases are smaller than 1K, so use a size that really is too small
edit('tests/preload_size_multi_db.sh', [
    ("run tiny --preload-size 1K", "run tiny --preload-size 8"),
    ("--preload-size 1K", "--preload-size 8") if False else ("for size in 600K 800K 2M 100M; do", "for size in 600K 800K 2M 100M; do"),
])
