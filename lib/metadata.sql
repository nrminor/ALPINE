WITH table_0 AS (
  SELECT
    *
  FROM
    read_csv('genbank_metadata.tsv')
)
SELECT
  *
FROM
  table_0
WHERE
  CHAR_LENGTH("Isolate Collection date") >= 10
  AND CHAR_LENGTH("Isolate Collection date") <> 0
  AND CHAR_LENGTH("Geographic location") <> 0

-- Generated by PRQL compiler version:0.11.1 (https://prql-lang.org)
