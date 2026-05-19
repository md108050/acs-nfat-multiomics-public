# src/01_baseline_table.R

source("src/00_setup.R")

# 出力先
out_dir <- "results/tables/baseline"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1. データの読み込み ----
df_cohort <- read_parquet("data/meta/cohort_191.parquet")

# 2. 変数の定義 ----
listvars <- c(
  "sex",
  "age",
  "BMI",
  "HbA1c",
  "smoking",
  "alcohol",
  "ACTH",
  "cortisol",
  "DHEAS",
  "dst_cortisol",
  "midnight_cortisol",
  "U_cortisol",
  "TC",
  "TG",
  "HDLC",
  "LDLC"
)

listcat  <- c("sex", "smoking", "alcohol")
listcont <- setdiff(listvars, listcat)

# 3. Table 1 の作成 ----
tab1 <- CreateTableOne(
  vars       = listvars,
  factorVars = listcat,
  strata     = "Diagnosis_2groups",
  data       = df_cohort,
  includeNA  = FALSE,
  addOverall = FALSE,
  testExact  = fisher.test,
  testNonNormal = kruskal.test
)

# 4. 出力 ----
res <- print(
  tab1,
  showAllLevels = TRUE,
  nonnormal  = listcont,
  smd        = FALSE,
  explain    = TRUE,
  catDigits  = 1,
  contDigits = 1,
  pDigits    = 3,
  test       = TRUE,
  format     = "fp"
)

write_csv(as.data.frame(res), file.path(out_dir, "Table1_2groups.csv"))
