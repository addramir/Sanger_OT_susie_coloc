from pyspark.sql import DataFrame, SparkSession
import pandas as pd
import pyspark.sql.functions as f
from pyspark.sql.types import *
from pyspark.sql.window import Window
import numpy as np
import os


max_ram_to_use="30g"

def get_spark_session(max_ram_to_use):
    return (
        SparkSession.builder
        .master('local[*]')
        .config("spark.driver.memory", max_ram_to_use)
        .appName('spark')
        .getOrCreate()
    )

spark = get_spark_session(max_ram_to_use)


v2d = spark.read.parquet("gs://genetics-portal-dev-data/22.09.1/outputs/v2d_credset")
v2d_99=v2d.filter(f.col("is99_credset")==True)


df=v2d_99.filter(f.col("study_id").rlike("FIN"))
FIN=df.toPandas()
FIN.to_csv("~/FINN_cs.csv",index=False)

df=v2d_99.filter(f.col("study_id").rlike("SAIGE"))
SAIGE=df.toPandas()
SAIGE.to_csv("~/SAIGE_cs.csv",index=False)