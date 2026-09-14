# boldigger3
dockerfile to create container with https://github.com/DominikBuchner/BOLDigger3 python package

To run this modified version of boldigger3:

docker run \
  -v /your/project/folder:/data \
  -w /data \
  joschlag/boldigger3:3.0.3 \
  python -m boldigger3 identify PATH_TO_FASTA PATH_TO_DATABASE --db 1 --mode 1


BOLD is constantly changing their platform.
Restrictions of POST submissions to the id engine
and currently the necessity for a login to download snapshots makes it necessary to update programs that uses their id engine.
