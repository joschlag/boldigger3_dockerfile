# boldigger3
dockerfile to create container with https://github.com/DominikBuchner/BOLDigger3 python package

To run this modified version of boldigger3:

docker run \
  -v /your/project/folder:/data \
  -w /data \
  joschlag/boldigger3:2.2.0 \
  python -m boldigger3 identify file_for_identification.fasta --db 1 --mode 1 --db-dir relative/path/for/database-file/boldigger3_db

This modified version deals with BOLDs constant changes of there platform.
The restriction of POST submissions to the id engine
and currently the necessity for a login to download snapshots.
The last point is not solved yet, therefore this version was changed to use an old database (duckdb) file. 
