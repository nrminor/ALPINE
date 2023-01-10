#!/usr/bin/env python3

import pandas as pd, sys, os

# splices gisaid data by dates and 
#
# tsv_path - path to gisaid data
# startdate - date to start splice at (inclusive) in the form of year-mm-dd, e.g. 2020-01-21
# enddate - date to stop splice at (inclusive) in the form of year-mm-dd, e.g. 2020-01-21
# geography - string key to search the location metadata, e.g. California
# throws exception if startdate is later than enddate
def splice_df_by_date_geo(tsv_path, startdate, enddate, geography):
    # Load tsv file and initialize variables
    df = pd.read_table(tsv_path, dtype=str)
    startdate = pd.to_datetime(startdate)
    enddate = pd.to_datetime(enddate)
    
    # error case: enddate is earlier than startdate
    if(enddate < startdate): raise ValueError("End date is before start date.")
    
    # Converting dates of the dataframe. 
    # please see for pd.to_datetime documentation: https://pandas.pydata.org/docs/reference/api/pandas.to_datetime.html
    # Parameter format = "%Y-%m-%d"
    # Parameter exact = True means that formats outside of that format are rejected.
    # Parameter errors = "coerce" means incorrectly formated dates or dates that cannot be parsed are changed to NaT.
    # When doing NaT >= startdate, the result is False (so these all rows with incorrect rows are automatically removed)
    dates = pd.to_datetime(df["Collection date"],errors = "coerce",format = "%Y-%m-%d", exact = True)
    
    # Create filter: dates are later than start, dates are earlier than end, and location contains the keyword geography
    filter_df = (dates >= startdate) & (dates <= enddate) & (df["Location"].str.contains(geography))
    
    # Creates local tsv file, with empty values represented as empty strings.
    df[filter_df].to_csv("temp_filtered_gisaid.tsv", sep="\t", na_rep='')
    
    return os.path.join(os.getcwd(),"temp_filtered_gisaid.tsv")

splice_df_by_date_geo(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])