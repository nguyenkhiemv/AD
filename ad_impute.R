#rerun pipelines
library("dplyr")
source( 'new.fimpute.R')
# value.vars = c("bmi", "GMV_TIV_out", "GDSTot", "GMV", "DFCorr", "AFIBRILL", "DISN", "nb1_score", "NITE")
value.vars = c("bmi", "GMV", "DFCorr", "nb1_score")
time.var = "age_rounded"

# Only use subjects with at least 2 observations
tbl = table(ad_dat$PIDN)
ok = names(tbl[tbl > 1.5])
ad_dat = ad_dat[ad_dat$PIDN %in% ok,]

# Select relevant columns
data <- dplyr::select( ad_dat, c( "PIDN", all_of(time.var), value.vars ) ) %>% data.frame( )
for (col in names(data))
  data[[col]] = as.numeric(as.character(data[[col]]))

id.var = "id"
data$id = as.numeric(as.factor(data$PIDN))

# Remove rows with all value.vars equal to NA
data = data[rowSums(is.na(data[,value.vars]))!=length(value.vars),]

pp = 2*length(value.vars)
lr = 0.1
niter = 1000
lambdas = c(0, 0.1, 1, 10)
res = cv.fimpute(data,
                 value.vars,
                 time.var,
                 id.var,
                 pp = pp,
                 lr = lr,
                 niter = niter,
                 lambdas = lambdas)
plot(lambdas, res$loss)
bestLambda = lambdas[which.min(res$loss)]


print(paste("Lambda selected in cross-validation:",bestLambda))
model = fimpute.fit(data,
                    value.vars,
                    time.var,
                    id.var,
                    pp = pp,
                    lr = lr,
                    niter = niter,
                    lambda = bestLambda)
