# Simple bash script that launches the make_hist.py script on apposite directories

# $1 and $2 should be local directories WITHOUT THE CHARACTER "/" AT THE END
# $3 is the matrix dimension
python make_hist.py $1 $3
python make_hist.py $2 $3
