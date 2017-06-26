To use jupyter notebooks:

## Terminal 1 (sherlock):

`cd /share/PI/rondror/docking_code/notebooks/`

1) Start a server

`sbatch ./start_server_notebook.sh`

2) Check that the server has started

`squeue -u <you>`

3) Look in server_host.out to see what node your server is running on

`cat server_host.out`


## Terminal 2 (local machine):

The following commands should be run in a Unix shell (windows users: download Cygwin).

`PORT=$((30000+RANDOM%29999))`

`ssh -t -L 8893:localhost:$PORT <you>@sherlock.stanford.edu ssh -L $PORT:localhost:8893 <node-name>`

An example of `<node-name>` is `gpu-28-1`.


## In your browser:

`http://localhost:8893/tree?`

If a token is requested, look in server_host.out.
