# Configuration for sync.sh. This file is committed. To override a value on one
# machine, export the variable before you run the script.

# The S3 prefix that holds data/. The bucket has versioning enabled and all four
# public access blocks on.
FLOWCYTO_S3_URI="s3://wguesdon-flow-cytometry/everything-flow-cytometry"

# The storage class used on push. INTELLIGENT_TIERING moves an object that is not
# read for 30 days into a cheaper tier, and it charges no retrieval fee.
FLOWCYTO_STORAGE_CLASS="INTELLIGENT_TIERING"
