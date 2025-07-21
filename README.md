# boldigger3
dockerfile to create container with https://github.com/DominikBuchner/BOLDigger3 python package

To run this modified version of boldigger3:

docker run \
  -v /your/project/folder:/data \
  -w /data \
  joschlag/boldigger3:2.0.1 \
  python -m boldigger3 identify file_for_identification.fasta --db 1 --mode 1 --db-dir relative/path/for/database-file/boldigger3_db
