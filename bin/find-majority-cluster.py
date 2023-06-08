#!/usr/bin/env python3

import sys
import os
import polars as pl

# Define the input files
if os.path.islink(sys.argv[1]):
    cluster_table = os.readlink(sys.argv[1])
else:
    cluster_table = sys.argv[1]

# read in cluster table as data frame
cluster_df = pl.read_csv(cluster_table, has_header = False, separator = "\t")

# define a function that writes the largest cluster number to disk
def find_largest_cluster(cluster_df: pl.DataFrame):

    if cluster_df.select(pl.count())[0,0] == 2:

        majority_cluster = cluster_df['column_2'].unique()[0]
        majority_centroid = cluster_df['column_9'].unique()[0]

    else:

        cluster_data = cluster_df.filter(pl.col("column_1") == "C")
        count = 0
        for row in cluster_data.rows():

            cluster_count = row[2]

            if cluster_count > count:
                count = cluster_count
                majority_cluster = row[1]
                majority_centroid = row[8]
            else:
                continue

    # Save the majority cluster to an environmental variable that will be accessed by nextflow
    os.environ["majority_cluster"] = str(majority_cluster)

    # Save the majority cemtroid as am environmental variable
    os.environ["majority_centroid"] = str(majority_centroid)

# define a function that writes the number of clusters to disk
def count_clusters(cluster_df: pl.DataFrame):

    cluster_data = cluster_df.filter(pl.col("column_1") == "C")
    cluster_count = len(cluster_data['column_2'].unique().to_list())

    # Save the majority cluster to a file
    os.environ["cluster_count"] = str(cluster_count)

# run the functions
find_largest_cluster(cluster_df) ; count_clusters(cluster_df)